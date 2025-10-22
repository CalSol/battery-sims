%% 21700 cell setup
cylindricalGeometry = batteryCylindricalGeometry(simscape.Value(0.070,"m"), simscape.Value(0.0105,"m"));

cylindricalCell = batteryCell(cylindricalGeometry);
cylindricalCell.CellModelOptions.BlockParameters.thermal_port = "model";
cylindricalCell.Mass = simscape.Value(0.068, "kg");
cylindricalCell.Capacity = simscape.Value(4.5, "A*hr");
cylindricalCell.Energy = simscape.Value(18.5, "W*hr");

%% 9p module (18 modules in series = 36s9p pack)
singleModule = batteryParallelAssembly(cylindricalCell, 9, 'ModelResolution', 'Detailed');

%big pack for thermal analysis, workaround for missing series assembly


packRepresentation = batteryParallelAssembly(cylindricalCell, ...
    18, ...  % 18 cells in parallel (2s9p modules x 9)
    'ModelResolution', 'Detailed');

fprintf('   Cells in representation: %d\n', packRepresentation.NumParallelCells);

%% Step 5: Display Pack Specifications
fprintf('\n=== 36s9p BATTERY PACK SPECIFICATIONS ===\n');
fprintf('Target Configuration: 36s9p (36 series × 9 parallel)\n');
fprintf('Module Arrangement: 18 modules, each 2s9p\n');
fprintf('Module Layout: 9 cells long × 2 cells wide per module\n');
fprintf('Cell type: 21700 Li-ion\n');
fprintf('Total cells: %d\n', 36 * 9);
fprintf('Pack nominal voltage: %.1fV (36 × 3.7V)\n', 36 * 3.7);
fprintf('Pack capacity: %.1fAh (9 × 5Ah)\n', 9 * 4.5);
fprintf('Pack energy: %.1fkWh\n', 36 * 9 * 18.5 / 1000);
fprintf('Pack mass: %.1fkg\n', 36 * 9 * 0.068);
fprintf('Energy density: %.0fWh/kg\n', (36 * 9 * 18.5) / (36 * 9 * 0.068));

fprintf('\nPer Module (2s9p):\n');
fprintf('Module voltage: %.1fV (2 × 3.7V)\n', 2 * 3.7);
fprintf('Module capacity: %.1fAh (9 × 5Ah)\n', 9 * 4.5);
fprintf('Module energy: %.1fWh (2 × 9 × 18.5Wh)\n', 2 * 9 * 18.5);
fprintf('Module mass: %.1fkg (18 × 68g)\n', 18 * 0.068);

%% top-down 36s9p cell grid simulation (no modules)
close all;
try
    % Grid size: 36 series by 9 parallel = 36 columns x 9 rows (top-down)
    nCols = 36; % series strings (x direction)
    nRows = 9;  % cells in parallel per string (y direction)

    % Thermal + electrical parameters (per cell)
    T_amb = 25;                % degC
    C_cell = 120;              % J/K thermal mass per cell
    R_internal = 0.012;        % Ohm per cell

    % Neighbor conductive couplings [W/K]
    Gx = 0.20;  % left/right conductance (between series neighbors)
    Gy = 0.20;  % up/down conductance (between parallel neighbors)

    % External cooling to ambient (lumped per cell)
    G_conv = 0.05;   % W/K (tunable convective/radiative loss)

    % Operating profile (per-cell C-rate)
    C_rate = 1.0;              % per-cell C rate
    I_cell = C_rate * 4.5;     % A per cell (4.5Ah)
    P_gen = I_cell^2 * R_internal; % W per cell, constant here

    % Time discretization
    t_grid = 0:10:3600;  % seconds, 10 s step for stability
    dt = t_grid(2) - t_grid(1);

    % State initialization (degC)
    T = T_amb * ones(nRows, nCols);

    % Diagnostics over time
    T_min = zeros(size(t_grid));
    T_max = zeros(size(t_grid));
    T_avg = zeros(size(t_grid));
    T_min(1) = min(T, [], 'all');
    T_max(1) = max(T, [], 'all');
    T_avg(1) = mean(T, 'all');

    % Explicit time stepping of 2D network: each cell exchanges with 4 neighbors + ambient loss
    for k = 2:numel(t_grid)
        Tprev = T;
        % Loop over grid (compact, readable; could be vectorized if needed)
        for r = 1:nRows
            for c = 1:nCols
                Tc = Tprev(r,c);
                % neighbors (use ambient when out of bounds)
                T_left  = (c > 1)     * Tprev(r, max(c-1,1))   + (c == 1)    * T_amb;
                T_right = (c < nCols) * Tprev(r, min(c+1,nCols)) + (c == nCols)* T_amb;
                T_up    = (r > 1)     * Tprev(max(r-1,1), c)   + (r == 1)    * T_amb;
                T_down  = (r < nRows) * Tprev(min(r+1,nRows), c) + (r == nRows)* T_amb;

                % Net conductive power into cell from neighbors
                P_cond = Gx*(T_left - Tc) + Gx*(T_right - Tc) + Gy*(T_up - Tc) + Gy*(T_down - Tc);
                % Lumped convection/radiation to ambient
                P_amb = G_conv*(T_amb - Tc);
                % Energy balance
                dT = (P_gen + P_cond + P_amb) * dt / C_cell;
                T(r,c) = Tc + dT;
            end
        end
        T_min(k) = min(T, [], 'all');
        T_max(k) = max(T, [], 'all');
        T_avg(k) = mean(T, 'all');
    end

    % Visualization: top-down heat map and simple diagnostics
    f = figure('Color','w','Position',[50 50 1700 600]);
    f.Name = '36s9p top-down thermal simulation';

    % 1) Pack summary traces
    subplot(1,3,1);
    plot(t_grid/60, T_avg, 'k', 'LineWidth', 2); hold on;
    plot(t_grid/60, T_min, '--', 'Color',[0.2 0.6 1.0], 'LineWidth', 1.5);
    plot(t_grid/60, T_max, '--', 'Color',[1.0 0.2 0.2], 'LineWidth', 1.5);
    yline(T_amb, ':', 'Ambient', 'Color',[0.3 0.3 0.3]);
    grid on;
    xlabel('Time (min)'); ylabel('Temperature (°C)');
    title(sprintf('36s9p pack: %.1fC per cell, R=%.3fΩ', C_rate, R_internal));
    legend('Average','Min','Max','Location','northwest');

    % 2) Final top-down heat map (rows=parallel, cols=series)
    subplot(1,3,2);
    imagesc(T); axis image; colormap(hot); colorbar; caxis([min(T_amb, min(T,[],"all")) max(T,[],"all")]);
    title(sprintf('Top-down temperature (final @ %.0f min)', t_grid(end)/60));
    xlabel('Series position (1→36)'); ylabel('Parallel position (1→9)');

    % 3) Column-wise (series string) average temperature bar plot
    subplot(1,3,3);
    colAvg = mean(T,1);
    bar(1:nCols, colAvg, 0.9, 'FaceColor', [0.2 0.6 1.0]); grid on;
    xlabel('Series string'); ylabel('Avg Temp (°C)');
    title('Series string average temperatures');

    % Package results for saving/analysis
    topdown = struct();
    topdown.params = struct('nRows', nRows, 'nCols', nCols, 'T_amb', T_amb, ...
                            'C_cell', C_cell, 'R_internal', R_internal, ...
                            'Gx', Gx, 'Gy', Gy, 'G_conv', G_conv, ...
                            'C_rate', C_rate, 'I_cell', I_cell, 'P_gen', P_gen);
    topdown.t = t_grid;
    topdown.T_final = T;
    topdown.T_min = T_min; topdown.T_max = T_max; topdown.T_avg = T_avg;

catch ME
    warning(ME.identifier, 'Top-down simulation failed: %s', ME.message);
end

%% simscape model to see circuit behavior
fprintf('\nBuilding Simscape models...\n');

% cleanup before build 
force_rebuild = true;  % set flag depending on your vibes
if force_rebuild
    pkgDirs = {'+singleGroup_9p', '+packThermal_18cell'};
    for k = 1:numel(pkgDirs)
        pkgDir = pkgDirs{k};
        if exist(pkgDir, 'dir')
            try
                rmdir(pkgDir, 's');
                fprintf('Deleted existing package directory: %s\n', pkgDir);
            catch ME
                warning('Could not delete package directory %s: %s', pkgDir, ME.message);
            end
        end
    end
   
    matArtifacts = {'singleGroup_9p.mat','packThermal_18cell.mat'};
    for k = 1:numel(matArtifacts)
        mf = matArtifacts{k};
        if exist(mf,'file')
            try
                delete(mf);
                fprintf('Deleted existing artifact: %s\n', mf);
            catch ME
                warning('Could not delete %s: %s', mf, ME.message);
            end
        end
    end
    libs = {'singleGroup_9p_lib','packThermal_18cell_lib'};
    for k = 1:numel(libs)
        lib = libs{k};
        if bdIsLoaded(lib)
            close_system(lib, 0); % close without saving
        end
        slxFile = [lib '.slx'];
        if exist(slxFile, 'file')
            try
                delete(slxFile);
                fprintf('Deleted existing library: %s\n', slxFile);
            catch ME
                warning('Could not delete %s: %s', slxFile, ME.message);
            end
        end
    end
end

if force_rebuild || ~exist('singleGroup_9p_lib.slx','file')
    buildBattery(singleModule,'LibraryName','singleGroup_9p');
    fprintf('\u2713 Single 9p group library built: singleGroup_9p_lib.slx\n');
end

% pack representation
if force_rebuild || ~exist('packThermal_18cell_lib.slx','file')
    buildBattery(packRepresentation,'LibraryName','packThermal_18cell');
    fprintf('\u2713 Pack representation library built: packThermal_18cell_lib.slx\n');
end

%% libraries for simscape
try
    if exist('singleGroup_9p_lib.slx', 'file')
        open_system('singleGroup_9p_lib');
    end
    
    if exist('packThermal_18cell_lib.slx', 'file')
        open_system('packThermal_18cell_lib');
    end
    
catch ME
    warning('library opening error - try closing matlab windows');
end

%% Step 9: simulate
fprintf('\n=== THERMAL ANALYSIS ===\n');

% thermal parameters
thermal_params = struct();
thermal_params.cell_thermal_mass = 120;  % J/K 
thermal_params.cell_thermal_resistance = 4.5;  % K/W 
thermal_params.internal_resistance = 0.012;  % Ohm 
thermal_params.ambient_temp = 25;  % °C
thermal_params.max_temp = 60;  % °C 

fprintf('Thermal Parameters:\n');
fprintf('  Cell thermal mass: %.0f J/K\n', thermal_params.cell_thermal_mass);
fprintf('  Cell thermal resistance: %.1f K/W\n', thermal_params.cell_thermal_resistance);
fprintf('  Internal resistance: %.3f Ohm\n', thermal_params.internal_resistance);
fprintf('  Ambient temperature: %.0f°C\n', thermal_params.ambient_temp);
fprintf('  Maximum safe temperature: %.0f°C\n', thermal_params.max_temp);


discharge_rates = [0.14, 0.25, 0.5, 1.0, 2.0, 3.0];  % C-rates
time_sim = 0:60:3600;  % 1 hour

figure('Name', 'Thermal Simulation Results', 'Position', [200 200 1200 800]);

for i = 1:length(discharge_rates)
    C_rate = discharge_rates(i);
    current = C_rate * 4.5; % 4.5ah capacity 
    power_per_cell = current^2 * thermal_params.internal_resistance; % W

    temp_rise = zeros(size(time_sim)); % rise above ambient (C)
    for t = 2:length(time_sim)
        dt = time_sim(t) - time_sim(t-1);
        heat_input = power_per_cell * dt;            % Joules
        temp_rise(t) = temp_rise(t-1) + heat_input / thermal_params.cell_thermal_mass; % no cooling term
    end

    temperature = thermal_params.ambient_temp + temp_rise;
    
    subplot(3,2,i);
    plot(time_sim/60, temperature, 'LineWidth', 2);
    hold on;
    plot([0 60], [thermal_params.max_temp thermal_params.max_temp], 'r--', 'LineWidth', 2);
    title(sprintf('%.1fC Discharge (%.1fA)', C_rate, current));
    xlabel('Time (minutes)');
    ylabel('Temperature (°C)');
    legend('Cell Temperature', 'Max Safe Temp', 'Location', 'southeast');
    grid on;
    
    
    max_temp_reached = max(temperature);
    if max_temp_reached > thermal_params.max_temp
        disp('too hot');
    end
end

% (Module-based conduction visualization removed to focus on cell-level top-down model)

% Save results (include top-down if available)
try
    if exist('topdown','var')
        save('thermal_simulation_results.mat', 'thermal_params', 'discharge_rates', 'time_sim', 'topdown');
    else
        save('thermal_simulation_results.mat', 'thermal_params', 'discharge_rates', 'time_sim');
    end
    fprintf('Results saved to: thermal_simulation_results.mat\n');
catch ME
    warning(ME.identifier, 'Could not save results: %s', ME.message);
end