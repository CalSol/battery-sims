

safe_mode = false;
show_plots = true;
enable_dynamic_visualization = ~safe_mode;

if safe_mode
    try
        s = settings;
        if hasGroup(s,'matlab') && hasGroup(s.matlab,'graphics') && hasSetting(s.matlab.graphics,'useSoftwareOpenGL')
            s.matlab.graphics.useSoftwareOpenGL.PersonalValue = true;
        end
    catch ME
        warning(ME.identifier, '%s', ME.message);
    end
end

%% Cell Definition (21700 Li-ion)
cylindricalGeometry = batteryCylindricalGeometry(simscape.Value(0.070,"m"), simscape.Value(0.0105,"m"));

cylindricalCell = batteryCell(cylindricalGeometry);
cylindricalCell.CellModelOptions.BlockParameters.thermal_port = "model";
cylindricalCell.Mass = simscape.Value(0.068, "kg");
cylindricalCell.Capacity = simscape.Value(4.5, "A*hr");
cylindricalCell.Energy = simscape.Value(18.5, "W*hr");

%% Pack Architecture (36s9p)
singleModule = batteryParallelAssembly(cylindricalCell, 9, 'ModelResolution', 'Detailed');
packRepresentation = batteryParallelAssembly(cylindricalCell, 18, 'ModelResolution', 'Detailed');

fprintf('\n36s9p Pack Specifications:\n');
fprintf('  Total cells: %d (36 series × 9 parallel)\n', 324);
fprintf('  Nominal voltage: %.1f V\n', 133.2);
fprintf('  Capacity: %.1f Ah\n', 40.5);
fprintf('  Energy: %.1f kWh\n', 5.40);
fprintf('  Mass: %.1f kg\n', 22.0);
fprintf('  Energy density: %.0f Wh/kg\n', 245);

%% Grid Thermal Simulation 
close all;
try
    nCols = 36;  % Series strings
    nRows = 9;   % Parallel cells per 
    T_amb = 25;
    
    % Cell geometry (21700)
    cell_diameter = 0.021;  % m (21mm)
    cell_height = 0.070;    % m (70mm)
    cell_mass = 0.068;      % kg
    
    % Thermal properties
    c_p = 1000;  % J/(kg·K) specific heat capacity for Li-ion
    C_cell = cell_mass * c_p;  % J/K thermal capacitance (~68 J/K)
    
    fprintf('\nCorrected Thermal Model Parameters:\n');
    fprintf('  Cell thermal capacitance: %.1f J/K\n', C_cell);
    
    contact_fraction = 0.08;  % 8% estimate
    contact_area = cell_height * cell_diameter * contact_fraction;  % m²
    
    % Cell-to-cell spacing
    cell_spacing = 0.002;  % 2mm gap
        
    k_contact = 2.0;   % W/(m·K)
    k_air = 0.026;     % W/(m·K)

    % Parallel thermal resistance model: contact + air gap
    A_contact = contact_area;
    A_air = (cell_height * cell_diameter) - contact_area;
    
    % Conductance = 1/R_thermal, where R = L/(k*A)
    G_contact = k_contact * A_contact / cell_spacing;
    G_air = k_air * A_air / cell_spacing;
    
    Gx = G_contact + G_air;  % W/K, series 
    Gy = Gx * 0.8;  % W/K, parallel 
    
    fprintf('  Cell-cell conductance (X): %.3f W/K\n', Gx);
    fprintf('  Cell-cell conductance (Y): %.3f W/K\n', Gy);
    
    % Convection to ambient
    h_conv = 8;  % W/(m²·K) 
    
    A_surface = pi * cell_diameter * cell_height;  % m²
    
    G_conv_interior = h_conv * A_surface * 0.3;  % Only 30% exposed for interior
    G_conv_edge = h_conv * A_surface * 0.5;      % 50% exposed for edges
    G_conv_corner = h_conv * A_surface * 0.7;    % 70% exposed for corners
    
    fprintf('  Convection conductance (interior): %.3f W/K\n', G_conv_interior);
    fprintf('  Convection conductance (edge): %.3f W/K\n', G_conv_edge);
    fprintf('  Convection conductance (corner): %.3f W/K\n', G_conv_corner);
    
    tau_thermal = C_cell / (2*Gx + 2*Gy + G_conv_interior);
    fprintf('  Thermal time constant (interior): %.0f s (%.1f min)\n', ...
            tau_thermal, tau_thermal/60);
    
    R_internal_base = 0.012;  
    C_rate = 1.0;
    I_cell = C_rate * 4.5;  % Amps
    
 
    T_ref = 25;  % °C reference temperature
    alpha_R = 0.005;  % 1/°C temperature coefficient
    
    % time
    t_end = 3600;  % seconds (1 hour)
    dt = 10;       % seconds time step
    t_grid = 0:dt:t_end;
    
    
    T = T_amb * ones(nRows, nCols);
    T_history = zeros(nRows, nCols, length(t_grid));
    T_history(:,:,1) = T;
    
    T_min = zeros(size(t_grid));
    T_max = zeros(size(t_grid));
    T_avg = zeros(size(t_grid));
    T_min(1) = T_amb;
    T_max(1) = T_amb;
    T_avg(1) = T_amb;
    
    fprintf('\nRunning thermal simulation...\n');
    fprintf('  Time step: %.1f s\n', dt);
    fprintf('  Total time: %.0f min\n', t_end/60);
    fprintf('  Discharge rate: %.1fC (%.1f A per cell)\n', C_rate, I_cell);
    
    % Time integration loop
    for k = 2:length(t_grid)
        Tprev = T;
        
        for r = 1:nRows
            for c = 1:nCols
                Tc = Tprev(r,c);
                
                % 1. Heat generation (temperature-dependent)
                R_cell = R_internal_base * (1 + alpha_R * (Tc - T_ref));
                P_gen = I_cell^2 * R_cell;  % Watts
                
                % 2. Conduction to neighboring cells
                if c > 1
                    T_left = Tprev(r, c-1);
                else
                    T_left = T_amb;  % Left boundary sees ambient
                end
                
                if c < nCols
                    T_right = Tprev(r, c+1);
                else
                    T_right = T_amb;  % Right boundary sees ambient
                end
                
                if r > 1
                    T_up = Tprev(r-1, c);
                else
                    T_up = T_amb;  % Top boundary sees ambient
                end
                
                if r < nRows
                    T_down = Tprev(r+1, c);
                else
                    T_down = T_amb;  % Bottom boundary sees ambient
                end
                
                Q_dot_left = Gx * (T_left - Tc);
                Q_dot_right = Gx * (T_right - Tc);
                Q_dot_up = Gy * (T_up - Tc);
                Q_dot_down = Gy * (T_down - Tc);
                
                P_cond = Q_dot_left + Q_dot_right + Q_dot_up + Q_dot_down;
                
                
                is_edge_r = (r == 1 || r == nRows);
                is_edge_c = (c == 1 || c == nCols);
                
                if is_edge_r && is_edge_c
                    G_conv_cell = G_conv_corner;  % Corner cell
                elseif is_edge_r || is_edge_c
                    G_conv_cell = G_conv_edge;    % Edge cell
                else
                    G_conv_cell = G_conv_interior; % Interior cell
                end
                
                P_conv = G_conv_cell * (T_amb - Tc);  % Watts (negative when Tc > T_amb)
                
                P_net = P_gen + P_cond + P_conv;
                dT = (P_net * dt) / C_cell;
                
                % Update temperature
                T(r,c) = Tc + dT;
            end
        end
        
        T_history(:,:,k) = T;
        T_min(k) = min(T(:));
        T_max(k) = max(T(:));
        T_avg(k) = mean(T(:));
        
        if mod(k, round(length(t_grid)/10)) == 0
            fprintf('  Progress: %.0f%% (t=%.0f min, T_avg=%.1f°C)\n', ...
                    100*k/length(t_grid), t_grid(k)/60, T_avg(k));
        end
    end
    
    fprintf('Simulation complete!\n');
    fprintf('  Final average temperature: %.1f°C\n', T_avg(end));
    fprintf('  Final temperature range: %.1f - %.1f°C\n', T_min(end), T_max(end));
    fprintf('  Maximum temperature gradient: %.1f°C\n', T_max(end) - T_min(end));
    
    % steady-state
    P_gen_avg = I_cell^2 * R_internal_base;
    G_total_avg = 2*Gx + 2*Gy + G_conv_interior;
    T_ss_predicted = T_amb + P_gen_avg / G_total_avg;
    fprintf('  Predicted steady-state (interior): %.1f°C\n', T_ss_predicted);
    
    %% Visualization
    figure('Color','w','Position',[50 50 1700 600]);
    
    % Temperature evolution
    subplot(1,3,1);
    plot(t_grid/60, T_avg, 'k', 'LineWidth', 2.5); hold on;
    plot(t_grid/60, T_min, '--', 'Color',[0.2 0.6 1.0], 'LineWidth', 1.5);
    plot(t_grid/60, T_max, '--', 'Color',[1.0 0.2 0.2], 'LineWidth', 1.5);
    yline(T_amb, ':', 'Ambient', 'Color',[0.3 0.3 0.3], 'LineWidth', 1.5);
    yline(60, '-.', 'Limit', 'Color',[0.8 0.2 0.2], 'LineWidth', 1.5);
    grid on; 
    xlabel('Time (min)', 'FontSize', 12); 
    ylabel('Temperature (°C)', 'FontSize', 12);
    title(sprintf('Temperature Evolution (%.1fC discharge)', C_rate), 'FontSize', 13);
    legend('Average','Min','Max','Location','southeast', 'FontSize', 10);
    xlim([0 t_end/60]);
    
    % Spatial temperature distribution
    subplot(1,3,2);
    imagesc(T); 
    axis image; 
    colormap(hot); 
    cb = colorbar;
    ylabel(cb, 'Temperature (°C)', 'FontSize', 11);
    clim([T_amb, max(T(:))]);
    title('Final Temperature Distribution', 'FontSize', 13);
    xlabel('Series String (1-36)', 'FontSize', 12); 
    ylabel('Parallel Position (1-9)', 'FontSize', 12);
    
    % Temperature by series string
    subplot(1,3,3);
    T_series_avg = mean(T, 1);  % Average across parallel strings
    T_series_min = min(T, [], 1);
    T_series_max = max(T, [], 1);
    
    plot(1:nCols, T_series_avg, 'k-o', 'LineWidth', 2, 'MarkerSize', 4); hold on;
    fill([1:nCols, nCols:-1:1], [T_series_min, fliplr(T_series_max)], ...
         [0.8 0.8 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    yline(T_amb, ':', 'Ambient', 'Color',[0.3 0.3 0.3]);
    grid on;
    xlabel('Series String Position', 'FontSize', 12);
    ylabel('Temperature (°C)', 'FontSize', 12);
    title('Series String Temperature Profile', 'FontSize', 13);
    legend('Average', 'Min-Max Range', 'Location', 'best', 'FontSize', 10);
    xlim([1 nCols]);
    
    % Store results
    topdown = struct();
    topdown.params = struct('nRows', nRows, 'nCols', nCols, ...
                            'T_amb', T_amb, 'C_cell', C_cell, ...
                            'R_internal_base', R_internal_base, ...
                            'alpha_R', alpha_R, ...
                            'Gx', Gx, 'Gy', Gy, ...
                            'G_conv_interior', G_conv_interior, ...
                            'G_conv_edge', G_conv_edge, ...
                            'G_conv_corner', G_conv_corner, ...
                            'C_rate', C_rate, 'I_cell', I_cell, ...
                            'tau_thermal', tau_thermal);
    topdown.t = t_grid;
    topdown.T_history = T_history;
    topdown.T_final = T;
    topdown.T_min = T_min;
    topdown.T_max = T_max;
    topdown.T_avg = T_avg;
    
catch ME
    warning(ME.identifier, 'Grid simulation failed: %s', ME.message);
    fprintf('Error location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% Simscape Model Generation
fprintf('\nBuilding Simscape models...\n');

force_rebuild = ~safe_mode;
if force_rebuild
    pkgDirs = {'+singleGroup_9p', '+packThermal_18cell'};
    for k = 1:numel(pkgDirs)
        if exist(pkgDirs{k}, 'dir')
            try
                rmdir(pkgDirs{k}, 's');
            catch ME
                warning('Could not delete %s: %s', pkgDirs{k}, ME.message);
            end
        end
    end
   
    matArtifacts = {'singleGroup_9p.mat','packThermal_18cell.mat'};
    for k = 1:numel(matArtifacts)
        if exist(matArtifacts{k},'file')
            try
                delete(matArtifacts{k});
            catch ME
                warning('Could not delete %s: %s', matArtifacts{k}, ME.message);
            end
        end
    end
    
    libs = {'singleGroup_9p_lib','packThermal_18cell_lib'};
    for k = 1:numel(libs)
        if bdIsLoaded(libs{k})
            close_system(libs{k}, 0);
        end
        slxFile = [libs{k} '.slx'];
        if exist(slxFile, 'file')
            try
                delete(slxFile);
            catch ME
                warning('Could not delete %s: %s', slxFile, ME.message);
            end
        end
    end
end

if force_rebuild || ~exist('singleGroup_9p_lib.slx','file')
    buildBattery(singleModule,'LibraryName','singleGroup_9p');
end

if force_rebuild || ~exist('packThermal_18cell_lib.slx','file')
    buildBattery(packRepresentation,'LibraryName','packThermal_18cell');
end

try
    if exist('singleGroup_9p_lib.slx', 'file')
        open_system('singleGroup_9p_lib');
    end
    if exist('packThermal_18cell_lib.slx', 'file')
        open_system('packThermal_18cell_lib');
    end
catch ME
    warning(ME.identifier, '%s', ME.message);
end

%% Thermal Analysis
fprintf('\nThermal Analysis\n');

thermal_params = struct();
thermal_params.cell_thermal_mass = 120;
thermal_params.cell_thermal_resistance = 4.5;
thermal_params.internal_resistance = 0.012;
thermal_params.ambient_temp = 25;
thermal_params.max_temp = 60;

discharge_rates = [0.14, 0.25, 0.5, 1.0, 2.0, 3.0];
time_sim = 0:60:3600;

if show_plots
    figure('Name', 'Thermal Simulation Results', 'Position', [200 200 1200 800]);
end

for i = 1:length(discharge_rates)
    C_rate = discharge_rates(i);
    current = C_rate * 4.5;
    power_per_cell = current^2 * thermal_params.internal_resistance;

    temp_rise = zeros(size(time_sim));
    for t = 2:length(time_sim)
        dt = time_sim(t) - time_sim(t-1);
        heat_input = power_per_cell * dt;
        heat_dissipated = temp_rise(t-1) / thermal_params.cell_thermal_resistance * dt;
        net_heat = heat_input - heat_dissipated;
        temp_rise(t) = temp_rise(t-1) + net_heat / thermal_params.cell_thermal_mass;
    end

    temperature = thermal_params.ambient_temp + temp_rise;
    
    if show_plots
        subplot(3,2,i);
        plot(time_sim/60, temperature, 'LineWidth', 2); hold on;
        plot([0 60], [thermal_params.max_temp thermal_params.max_temp], 'r--', 'LineWidth', 2);
        title(sprintf('%.1fC (%.1fA)', C_rate, current));
        xlabel('Time (min)'); ylabel('Temperature (°C)');
        legend('Cell Temp', 'Max Safe', 'Location', 'southeast');
        grid on;
    end
    
    max_temp = max(temperature);
    steady_state = thermal_params.ambient_temp + (power_per_cell * thermal_params.cell_thermal_resistance);
    
    fprintf('  %.2fC: Max=%.1f°C, SS=%.1f°C', C_rate, max_temp, steady_state);
    if max_temp > thermal_params.max_temp
        fprintf(' [EXCEEDS LIMIT]\n');
    else
        fprintf('\n');
    end
end

try
    if exist('topdown','var')
        save('thermal_simulation_results.mat', 'thermal_params', 'discharge_rates', 'time_sim', 'topdown');
    else
        save('thermal_simulation_results.mat', 'thermal_params', 'discharge_rates', 'time_sim');
    end
catch ME
    warning(ME.identifier, 'Could not save results: %s', ME.message);
end

%% Dynamic Visualization with Custom Thermal Model Integration
dynviz_should_run = enable_dynamic_visualization;
try
    if evalin('base','exist(''dynviz_inprogress'',''var'')') && evalin('base','dynviz_inprogress')
        if bdIsLoaded('packVisualizationModel')
            dynviz_should_run = false;
        else
            evalin('base','clear dynviz_inprogress');
        end
    end
catch
end

if dynviz_should_run
    assignin('base','dynviz_inprogress',true);
    dynviz_cleanup = onCleanup(@() assignin('base','dynviz_inprogress',false)); %#ok<NASGU>
    try
        fprintf('\nDynamic Visualization with Custom Thermal Integration\n');
        
        % Use actual cell geometry from earlier definition
        vizCylGeom = batteryCylindricalGeometry(simscape.Value(0.070,"m"), simscape.Value(0.0105,"m"));
        vizCell = batteryCell(vizCylGeom);
        vizCell.CellModelOptions.CellModelBlockPath = "batt_lib/Cells/Battery Equivalent Circuit";
        vizCell.CellModelOptions.BlockParameters.ThermalModel = "LumpedThermalMass";
        vizCell.Mass = simscape.Value(0.068, "kg");
        vizCell.Capacity = simscape.Value(4.5, "A*hr");
        
        % Create pack structure matching 36s9p layout (9 parallel, 36 series)
        vizPA = batteryParallelAssembly(vizCell, 9, ModelResolution="Detailed", Topology="Hexagonal");
        vizModule = batteryModule(vizPA, 36, StackingAxis="X", ...
                                   InterParallelAssemblyGap=simscape.Value(0.005,"m"), ...
                                   ModelResolution="Detailed");
        vizMA = batteryModuleAssembly(vizModule, InterModuleGap=simscape.Value(0.001,"m"), StackingAxis="X");
        pack = batteryPack(vizMA, AmbientThermalPath="CellBasedThermalResistance", ...
                                  CoolantThermalPath="CellBasedThermalResistance");

        if show_plots
            f = figure(Color="w",Position=[100 100 1000 500]); 
            tl = tiledlayout(f,1,2,TileSpacing="Compact");
            nexttile; batteryChart(f, pack); title('36s9p Pack Layout');
            nexttile; batteryChart(f, pack, SimulationStrategyVisible="On"); title('Simulation Strategy');
        end

        libraryname = "packVisualizationDyn";
        
        if bdIsLoaded(libraryname + "_lib"), close_system(libraryname + "_lib",0); end
        if exist(libraryname + "_lib.slx",'file'), delete(libraryname + "_lib.slx"); end
        if exist(libraryname + ".slx",'file'), delete(libraryname + ".slx"); end
        if exist(libraryname + ".mat",'file'), delete(libraryname + ".mat"); end
        if exist(libraryname + "_param.m",'file'), delete(libraryname + "_param.m"); end
        pkgDir = "+" + libraryname;
        if exist(pkgDir, 'dir')
            try
                rmdir(pkgDir, 's');
            catch ME
                warning('Could not delete package directory %s: %s', pkgDir, ME.message);
            end
        end
        
        fprintf('Building battery pack library (36s9p = %d cells)...\n', nCols*nRows);
        buildBattery(pack, LibraryName=libraryname, MaskParameters="VariableNamesByType", MaskInitialTargets="VariableNamesByInstance");

        modelname = "packVisualizationModel";
        if bdIsLoaded(modelname), close_system(modelname,0); end
        if exist(modelname + ".slx",'file'), delete(modelname + ".slx"); end
        open_system(new_system(modelname));

        batteryBlockPath = strcat(modelname,"/",pack.Name);
        electricalRefBlockPath = strcat(modelname,"/","ElectricalReference");
        solverConfigBlockPath = strcat(modelname,"/","Solver");
        currentSourceBlockPath = strcat(modelname,"/","DC Current Source");
        ambientTemperatureSourceBlockPath = strcat(modelname,"/","Temperature Source");
        coolantTemperatureSourceBlockPath = strcat(modelname,"/","Temperature Source Clnt");

        add_block(strcat(libraryname,"/",pack.Name), batteryBlockPath, position=[150,100,300,250]);
        add_block("fl_lib/Electrical/Electrical Elements/Electrical Reference", electricalRefBlockPath, position=[165,320,185,340], orientation="down", ShowName="off");
        add_block("nesl_utility/Solver Configuration", solverConfigBlockPath, position=[-80,280,-30,320]);
        add_block("fl_lib/Thermal/Thermal Sources/Temperature Source", ambientTemperatureSourceBlockPath, position=[205,300,245,340], ShowName="off");
        add_block("fl_lib/Thermal/Thermal Sources/Temperature Source", coolantTemperatureSourceBlockPath, position=[255,300,295,340], ShowName="off");
        add_block("fl_lib/Electrical/Electrical Sources/DC Current Source", currentSourceBlockPath, position=[-30,150,0,200], orientation="down", i0=num2str(-81));

        batteryBlockPortHandles = get_param(batteryBlockPath,"PortHandles");
        electricalRefBlockPortHandles = get_param(electricalRefBlockPath,"PortHandles");
        solverConfigBlockPortHandles = get_param(solverConfigBlockPath,"PortHandles");
        currentSourceBlockPortHandles = get_param(currentSourceBlockPath,"PortHandles");
        ambientTemperatureSourceBlockPortHandles = get_param(ambientTemperatureSourceBlockPath,"PortHandles");
        coolantTemperatureSourceBlockPortHandles = get_param(coolantTemperatureSourceBlockPath,"PortHandles");

        add_line(modelname, batteryBlockPortHandles.RConn(1), currentSourceBlockPortHandles.RConn, autorouting="smart");
        add_line(modelname, batteryBlockPortHandles.LConn, currentSourceBlockPortHandles.LConn, autorouting="smart");
        add_line(modelname, batteryBlockPortHandles.RConn(1), electricalRefBlockPortHandles.LConn, autorouting="smart");
        add_line(modelname, batteryBlockPortHandles.RConn(1), solverConfigBlockPortHandles.RConn, autorouting="smart");
        add_line(modelname, batteryBlockPortHandles.RConn(2), ambientTemperatureSourceBlockPortHandles.LConn, autorouting="smart");
        add_line(modelname, batteryBlockPortHandles.RConn(3), coolantTemperatureSourceBlockPortHandles.LConn, autorouting="smart");

        set_param(modelname,"SimscapeLogType","all");
        run(strcat(libraryname,"_param.m"));

        % Map custom thermal grid results to Simscape pack structure
        fprintf('Mapping custom thermal grid to Simscape pack...\n');
        if exist('ModuleAssembly1','var') && isfield(ModuleAssembly1,'Module1')
            numCells = nRows * nCols;  % 324 cells
            
            T_vector_C = reshape(T', numCells, 1); 
            T_vector_K = T_vector_C + 273.15;
            
            fprintf('  Temperature range: %.1f°C to %.1f°C\n', min(T_vector_C), max(T_vector_C));
            fprintf('  Setting %d cell temperatures from grid simulation\n', numCells);
            
            % Initialize all cells with grid-based temperatures
            ModuleAssembly1.Module1.socCell = 0.8 * ones(numCells, 1);
            ModuleAssembly1.Module1.batteryTemperature = T_vector_K;
            
            % Set thermal resistance based on grid model
            % Convert conductance (W/K) to resistance (K/W): R = 1/G
            try
                ModuleType1 = evalin('base', 'ModuleType1');
                
                R_amb_avg = 1 / (2*Gx + 2*Gy + G_conv);
                ModuleType1.AmbientResistance = R_amb_avg * ones(1, numCells);
                
                % Coolant resistance (set high since we model ambient cooling)
                ModuleType1.CoolantResistance = 1e6 * ones(1, numCells);
                
                % Apply position-dependent thermal resistance
                for idx = 1:numCells
                    % Get row and column in grid
                    r = mod(idx-1, nRows) + 1;  % Parallel position (1-9)
                    c = floor((idx-1) / nRows) + 1;  % Series position (1-36)
                    
                    % Edge cells have better cooling (more exposed surfaces)
                    num_exposed = 0;
                    if r == 1 || r == nRows, num_exposed = num_exposed + 1; end
                    if c == 1 || c == nCols, num_exposed = num_exposed + 1; end
                    
                    % Reduce thermal resistance for edge cells
                    if num_exposed > 0
                        edge_factor = 0.7;  % 30% better cooling
                        ModuleType1.AmbientResistance(idx) = R_amb_avg * edge_factor^num_exposed;
                    end
                end
                
                fprintf('  Thermal resistance: %.1f to %.1f K/W\n', ...
                    min(ModuleType1.AmbientResistance), max(ModuleType1.AmbientResistance));
                
                assignin('base', 'ModuleType1', ModuleType1);
            catch ME
                warning(ME.identifier, 'Could not set thermal resistances: %s', ME.message);
            end
        end

        fprintf('Running Simscape simulation...\n');
        out = sim(modelname, StartTime="0", StopTime="600");
        
        % Create battery simulation log for visualization
        fprintf('Creating dynamic visualizations...\n');
        batterySimLog = batterySimulationLog(pack, out.simlog.(pack.Name));

        if show_plots
            % Temperature visualization
            f1 = uifigure('Color','w','Name','Cell Temperature Distribution');
            f1.Position = [50 400 1200 600];
            g1 = uigridlayout(f1,[1 1]);
            packChart = batterySimulationChart(g1, batterySimLog);
            packChartColorBar = colorbar(packChart);
            ylabel(packChartColorBar, strcat("Cell temperature"," (",batterySimLog.SelectedVariableUnit,")"), 'FontSize',14);

            % Current distribution visualization
            batteryCurrentSimLog = batterySimLog;
            batteryCurrentSimLog.SelectedVariable = "batteryCurrent";
            
            % State of charge visualization
            batteryStateOfChargeSimLog = batterySimLog;
            batteryStateOfChargeSimLog.SelectedVariable = "socCell";

            f2 = uifigure('Color','w','Name','Current and State of Charge');
            f2.Position = [50 50 1400 600];
            g2 = uigridlayout(f2,[1 2]);
            
            packCurrentChart = batterySimulationChart(g2, batteryCurrentSimLog);
            packCurrentChartColorBar = colorbar(packCurrentChart);
            ylabel(packCurrentChartColorBar, strcat("Cell current"," (",batteryCurrentSimLog.SelectedVariableUnit,")"), 'FontSize',14);
            
            packStateOfChargeChart = batterySimulationChart(g2, batteryStateOfChargeSimLog);
            packStateOfChargeChartColorBar = colorbar(packStateOfChargeChart);
            ylabel(packStateOfChargeChartColorBar, strcat("Cell state of charge"," (",batteryStateOfChargeSimLog.SelectedVariableUnit,")"), 'FontSize',14);
            
            % Compare custom thermal model with Simscape results
            figure('Name','Thermal Model Comparison','Position',[1300 400 600 600],'Color','w');
            
            subplot(2,1,1);
            imagesc(T); axis image; colormap(hot); colorbar;
            clim([min(T_amb, min(T,[],"all")) max(T,[],"all")]);
            title('Custom Grid Model - Final Temperature (°C)');
            xlabel('Series position (1-36)'); ylabel('Parallel position (1-9)');
            
            subplot(2,1,2);
            % Extract final temperatures from simulation
            T_simscape = ModuleAssembly1.Module1.batteryTemperature - 273.15;
            T_simscape_grid = reshape(T_simscape, [nRows, nCols])';
            imagesc(T_simscape_grid); axis image; colormap(hot); colorbar;
            clim([min(T_amb, min(T,[],"all")) max(T,[],"all")]);
            title('Simscape Initial Conditions from Grid (°C)');
            xlabel('Series position (1-36)'); ylabel('Parallel position (1-9)');
        end
        
        fprintf('\nVisualization complete!\n');
        fprintf('  Grid model final avg temp: %.1f°C\n', mean(T, 'all'));
        fprintf('  Simscape model cells: %d\n', pack.NumModels);
        fprintf('  Use the interactive charts to explore temperature evolution\n');
        
    catch ME
        warning(ME.identifier,'Visualization error: %s',ME.message);
        if exist('modelname','var') && bdIsLoaded(modelname)
            try
                close_system(modelname, 0);
            catch
            end
        end
    end
end