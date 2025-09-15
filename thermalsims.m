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

%% ui stuff
close all;
try
    f = figure('Color','w','Position',[50 50 1800 600]);
    f.Name = 'thermal analysis';
    
    % simulate temperature distribution across pack
    ambient_temp = 25; % °C
    num_modules = 18;
    num_cells_per_module = 9;
    discharge_rate = 1.0; % 1c
    current = discharge_rate * 4.5; % 4.5ah capacity
    power_per_cell = current^2 * 0.012; % 0.012 ir
    cell_thermal_mass = 120;
    cell_thermal_resistance = 4.5; % kept for reference, no longer used for cooling
    % --- radiation parameters (simplified) ---
    sigma_SB = 5.670374419e-8;  % boltzmann constant (W/m^2/K^4)
    emissivity = 0.8;           % assumed emissivity
    A_face = 0.003;             % area between neighboring cells? need to cad this out and verify

    time_sim = 0:60:3600;
    % store absolute temperatures in Kelvin for radiation calc
    temp_cells_K = (ambient_temp) * ones(num_cells_per_module, length(time_sim)) + 273.15; % initialize at ambient

    for t = 2:length(time_sim)
        dt = time_sim(t) - time_sim(t-1);
        heat_input = power_per_cell * dt;         
        for c = 1:num_cells_per_module
            T_prev = temp_cells_K(c, t-1);
            net_Q_rad_W = 0; % net power (W)
            if c > 1
                T_neighbor = temp_cells_K(c-1, t-1);
                net_Q_rad_W = net_Q_rad_W + emissivity * sigma_SB * A_face * (T_prev^4 - T_neighbor^4);
            end
            if c < num_cells_per_module
                T_neighbor = temp_cells_K(c+1, t-1);
                net_Q_rad_W = net_Q_rad_W + emissivity * sigma_SB * A_face * (T_prev^4 - T_neighbor^4);
            end
            % energy balance (thank you radke)
            dT = (heat_input - net_Q_rad_W * dt) / cell_thermal_mass;
            temp_cells_K(c, t) = T_prev + dT;
        end
    end
    % convert back to Celsius for plotting
    cell_temperatures = temp_cells_K - 273.15;
    % max cell temp in module taken to represent (have to fix this later)

    % Module temperature placeholders (will be updated after full thermal simulation later)
    module_temps = ambient_temp * ones(1, num_modules); % start at ambient
    temp_distribution = module_temps; % alias
    % We will update these rectangles after computing lumped thermal profiles (see later section).
    
    
    subplot(1,3,1);
    
    hold on;
    
    % Draw 18 modules in a single horizontal line
    cmap = jet(256); % color map
    tmin = min(temp_distribution);
    tmax = max(temp_distribution);
    rectHandles = gobjects(1, num_modules);
    tempTextHandles = gobjects(1, num_modules);
    for module = 1:num_modules
        x_pos = module * 2;
        y_pos = 5;
        % Map temperature to color (handle uniform case to avoid divide-by-zero)
        if tmax == tmin
            temp_normalized = 128; % mid-point color when all temps identical
        else
            temp_normalized = round(1 + (temp_distribution(module) - tmin) / (tmax - tmin) * 255);
            temp_normalized = min(max(temp_normalized,1),256); % Clamp to [1,256]
        end
        module_color = cmap(temp_normalized, :);
        % modules to see temperature gradient
        rectHandles(module) = rectangle('Position', [x_pos-0.8, y_pos-1.5, 1.6, 3], ...
                 'FaceColor', module_color, 'EdgeColor', 'black', 'LineWidth', 1.5);
        % label
        text(x_pos, y_pos, sprintf('M%d', module), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        % label temp
        tempTextHandles(module) = text(x_pos, y_pos-2.5, sprintf('%.1f°C', temp_distribution(module)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', 'red');
    end
    
    title('Pack Configuration: 18 Modules — Awaiting Simulation Update', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Module Position Along Pack Length');
    ylabel('Pack Width');
    xlim([0, 38]); ylim([0, 8]);
    axis equal;
    grid on;
    
    subplot(1,3,2);
    % single cell model -- making things much more advanced
    hold on;
    
    % single cell thermal visualization
    [theta, r] = meshgrid(linspace(0, 2*pi, 20), linspace(0, 1, 10));
    [x_cell, y_cell] = pol2cart(theta, r);
    
    % radial temperature distribution INSIDE cell
    temp_radial = ambient_temp + 15 * (1 - r.^2); 
    
    % plot distribution based on what math was above this line

    surf(x_cell, y_cell, zeros(size(x_cell)), temp_radial, 'EdgeColor', 'none');
    % thermal zones
    contour(x_cell, y_cell, temp_radial, [30, 35, 40], 'LineColor', 'black', 'LineWidth', 2);
    title('Single Cell Thermal Strategy', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Cell Radial Position');
    ylabel('Cell Radial Position');
    view(0, 90); % look at it
    axis equal;
    cb = colorbar;
    ylabel(cb, 'Temperature (°C)');
    
    % annotations that ai made
    text(0, -1.5, 'Core: 40°C', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
    text(0.7, 0.7, 'Edge: 30°C', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'blue');
    text(0, 1.5, '21700 Li-ion Cell\nRadial Heat Distribution', 'HorizontalAlignment', 'center', 'FontSize', 9);
    
    subplot(1,3,3);
    % PANEL 3: pack heatmap based on cell temps from above
    hold on;
    
    % temp zones based on module temps
    temp_zones = zeros(1, num_modules);
    for i = 1:num_modules
        if temp_distribution(i) <= 30
            temp_zones(i) = 1; % Cool - Blue
        elseif temp_distribution(i) <= 35
            temp_zones(i) = 2; % Medium - Yellow
        else
            temp_zones(i) = 3; % Hot - Red
        end
    end
    
    % colors that ARE BETWEEN [0, 1]
    zone_colors = [0.2 0.6 1.0;   % blue
                   1.0 0.9 0.2;   % yellow
                   1.0 0.2 0.2];  % red
    
    % heatmap
    for module = 1:18
        x_pos = module * 2;
        y_pos = 5;
        
        color_idx = temp_zones(module);
        
        
        rectangle('Position', [x_pos-0.8, y_pos-1.5, 1.6, 3], ...
                 'FaceColor', zone_colors(color_idx, :), 'EdgeColor', 'black', 'LineWidth', 2);
        
        text(x_pos, y_pos, sprintf('M%d', module), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        
        zone_labels = {'COOL', 'MED', 'HOT'};
        text(x_pos, y_pos-2.5, sprintf('%s\n%.1f°C', zone_labels{color_idx}, temp_distribution(module)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
    end
    
    title('Pack Thermal Heatmap: 3-Zone Strategy', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Module Position Along Pack Length');
    ylabel('Pack Width');
    xlim([0, 38]); ylim([0, 8]);
    axis equal;
    grid on;
    

    legend_x = 32;
    legend_y = 7;
    
    % Cool zone
    rectangle('Position', [legend_x, legend_y, 1, 0.5], 'FaceColor', zone_colors(1,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y + 0.25, 'COOL (≤30°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
    
    % Medium zone  
    rectangle('Position', [legend_x, legend_y-0.8, 1, 0.5], 'FaceColor', zone_colors(2,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y - 0.55, 'MEDIUM (31-35°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
    
    % Hot zone
    rectangle('Position', [legend_x, legend_y-1.6, 1, 0.5], 'FaceColor', zone_colors(3,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y - 1.35, 'HOT (>35°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
    

catch ME
    fprintf('you messed up something in the ui bro\n');
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
        print('too hot')
    end
end

% --- Per-module conduction model & heatmap update (1C case) ---
try
    target_C = 1.0;          % choose C-rate to visualize spatial distribution
    I1 = target_C * 4.5;     % A
    P_cell = I1^2 * thermal_params.internal_resistance; % W per cell
    cells_per_module = num_cells_per_module; % using earlier definition (9)
    P_module = cells_per_module * P_cell;    % uniform generation per module
    C_module = cells_per_module * thermal_params.cell_thermal_mass; % J/K per module

    % Conduction only (no external cooling): set inter-module conductance, remove ambient losses.
    G_cond = 0.18;   % W/K between adjacent modules (effective lumped conductance)

    t_vec = time_sim; % reuse 0:60:3600
    moduleTemps = thermal_params.ambient_temp * ones(length(t_vec), num_modules);
    % 1D conduction with fixed ambient boundary (ghost nodes at ambient).
    % This produces a spatial gradient: edges lose heat to ambient, center runs hotter.
    for ti = 2:length(t_vec)
        dt = t_vec(ti) - t_vec(ti-1);
        Tprev = moduleTemps(ti-1, :);
        Tnew = Tprev;
        for m = 1:num_modules
            Tm = Tprev(m);
            T_left  = thermal_params.ambient_temp; if m > 1, T_left  = Tprev(m-1); end
            T_right = thermal_params.ambient_temp; if m < num_modules, T_right = Tprev(m+1); end
            % Net conductive power INTO module m from neighbors & ambient-boundary nodes
            P_cond = G_cond * ((T_left - Tm) + (T_right - Tm)); % W
            net_power = P_module + P_cond; % generation plus conductive in-flow (can be negative)
            dT = (net_power * dt) / C_module;
            Tnew(m) = Tm + dT;
        end
        moduleTemps(ti,:) = Tnew;
    end

    finalTemps = moduleTemps(end, :);
    tmin = min(finalTemps); tmax = max(finalTemps);
    cmap = jet(256);
    for m = 1:num_modules
        if tmax == tmin
            ci = 128;
        else
            ci = round(1 + (finalTemps(m) - tmin)/(tmax - tmin)*255);
            ci = min(max(ci,1),256);
        end
        set(rectHandles(m), 'FaceColor', cmap(ci,:));
        set(tempTextHandles(m), 'String', sprintf('%.1f°C', finalTemps(m)));
    end
    subplot(1,3,1);
    title(sprintf('Per-Module Conduction (Edges Fixed at Ambient)  Min %.1f°C  Max %.1f°C', tmin, tmax), 'FontSize', 14, 'FontWeight', 'bold');
    drawnow;

    % --- Redraw zone heatmap (subplot 3) with new distribution ---
    subplot(1,3,3); cla; hold on;
    updatedTemps = finalTemps;
    temp_zones = zeros(1, num_modules);
    for iZ = 1:num_modules
        if updatedTemps(iZ) <= 30
            temp_zones(iZ) = 1;
        elseif updatedTemps(iZ) <= 35
            temp_zones(iZ) = 2;
        else
            temp_zones(iZ) = 3;
        end
    end
    zone_colors = [0.2 0.6 1.0; 1.0 0.9 0.2; 1.0 0.2 0.2];
    for module = 1:num_modules
        x_pos = module * 2; y_pos = 5; color_idx = temp_zones(module);
        rectangle('Position', [x_pos-0.8, y_pos-1.5, 1.6, 3], 'FaceColor', zone_colors(color_idx,:), 'EdgeColor', 'black', 'LineWidth', 2);
        text(x_pos, y_pos, sprintf('M%d', module), 'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        zone_labels = {'COOL','MED','HOT'};
        text(x_pos, y_pos-2.5, sprintf('%s\n%.1f°C', zone_labels{color_idx}, updatedTemps(module)), 'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
    end
    title(sprintf('Pack Thermal Zones (1C, Edge Ambient BC)  Max %.1f°C', max(updatedTemps)), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Module Position Along Pack Length'); ylabel('Pack Width'); xlim([0, 38]); ylim([0,8]); axis equal; grid on;

    % Legend
    legend_x = 32; legend_y = 7;
    rectangle('Position', [legend_x, legend_y, 1, 0.5], 'FaceColor', zone_colors(1,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y + 0.25, 'COOL (≤30°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
    rectangle('Position', [legend_x, legend_y-0.8, 1, 0.5], 'FaceColor', zone_colors(2,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y - 0.55, 'MED (31–35°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
    rectangle('Position', [legend_x, legend_y-1.6, 1, 0.5], 'FaceColor', zone_colors(3,:), 'EdgeColor', 'black');
    text(legend_x + 1.2, legend_y - 1.35, 'HOT (>35°C)', 'FontSize', 9, 'VerticalAlignment', 'middle');
catch ME
    warning(ME.identifier, '%s', ME.message);
end

save('thermal_simulation_results.mat', 'thermal_params', 'discharge_rates', 'time_sim');
fprintf('Results saved to: thermal_simulation_results.mat\n');