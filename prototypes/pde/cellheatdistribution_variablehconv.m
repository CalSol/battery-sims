% Cell geometry (m, kg).
D=0.02155;
r=D/2;
h=0.070155;
m=0.07;
V=pi*r^2*h;
SA=2*pi*(r^2)+2*pi*(2*r)*h;
rho = m/V;   % density, kg/m^3
% Electrical parameters.
V_nominal=3.6;
R_i=0.015;
% Discharge test conditions.
I=10;    %choose a discharge current, A
T_ambient=25+273.15;
T_0=T_ambient;
% Thermal properties.
C_p=1040;   % specific heat, J/(kg*K)
k_r=0.9;  % radial thermal conductivity, W/(m*K)
k_theta=k_r;    % cylindrical symmetry
k_z=24.2;    % axial thermal conductivity, W/(m*K)
q=(I^2)*(R_i)/V; % volumetric heat generation, W/m^3

% Transient heat equation: rho*C_p*dT/dt = div(k*grad(T)) + q.

thermalmodel = createpde('thermal','transient');
thermalmodel.Geometry = multicylinder(r, h);
figure; pdegplot(thermalmodel);
allfaces=1:thermalmodel.Geometry.NumFaces;
hmax=0.0025;
generateMesh(thermalmodel, 'Hmax',hmax);
figure; pdemesh(thermalmodel)
k=[k_r;k_theta;k_z];
% Define thermal properties
thermalProperties(thermalmodel, 'ThermalConductivity',k,'MassDensity',q,'SpecificHeat',C_p);
internalHeatSource(thermalmodel, q);
% Apply the spatially varying convection boundary condition and initial state.
thermalBC(thermalmodel, "Face",allfaces,"HeatFlux",@(region,state) h_flux(region, state, T_ambient, h_min, h_max));
thermalIC(thermalmodel,T_0);    

% h_conv as a ftn of radial distance from x-axis (assuming airflow coming
% in from x-axis)
h_min = 10;
h_max = 30;

function q_flux = h_flux(region, state, T_ambient, h_min, h_max)
    T_local = state.u.';
    theta = atan2(region.x, region.y);
    h_conv = h_min + (h_max - h_min).*((1+cos(theta))/2);   %.* means element-wise multiplication
    q_flux = (h_conv * (T_local - T_ambient)).';
end

% Solve the transient thermal model
tlist = 0:5:3600; % time vector
result = solve(thermalmodel, tlist);
T_f=result.Temperature(:,end);

%%
% Visualize the final temperature field.
figure;
pdeplot3D(thermalmodel,"ColorMapData",T_f-273.15)
title(sprintf("Heat distribution at I=10 A after 3600s, 25°C ambient T"))
colorbar;
%%
% Plot the radial temperature history at mid-height.
d_x=linspace(0,r,100);
d_y=zeros(size(linspace(0,r,100)));
d_z=ones(size(linspace(0,r,100)))*(h/2);

T_radial=zeros(100, length(tlist));
for i=1:length(tlist)
    T_radial(:,i) = interpolateTemperature(result, d_x, d_y, d_z, i);
end

figure;
surf(tlist,d_x,T_radial)
xlabel(sprintf("Time (s)"))
ylabel(sprintf("Distance from center axis (m)"))
zlabel(sprintf("Temperature (K)"))
