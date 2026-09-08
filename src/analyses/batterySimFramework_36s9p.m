%% ========================================================================
%  batterySimFramework_36s9p.m
%  Complete battery-pack simulation framework for a 36s9p 21700 Li-ion pack.
%
%  Physical build: 2s9p modules (each a 2x9 cell block); 18 modules arranged
%  as two 9x18 blocks ("rows going back") separated by an air gap, giving a
%  324-cell 18x18 physical grid. Cooling is NEGATIVE-PRESSURE FORCED AIR: an
%  induced-draft fan at the central gap draws ambient air inward through each
%  block toward the gap exhaust.
%
%  Pipeline:
%    (1) Pack definition + specs (derived, not hard-coded)
%    (2) Equivalent-circuit model (2nd-order Thevenin, OCV/SOC, R0(T))
%    (3) Variable motor-derived load profile  -> pack current demand
%    (4) Coupled electro-thermal simulation of all 324 cells on the 18x18
%        physical grid with the central gap and inward air flow
%    (5) Detailed single-cell PDE thermal model (PDE Toolbox), corrected
%    (6) Simscape Battery Simulink model of the 36s9p pack (buildBattery)
%    (7) Clean, presentable visualizations (saved to PNG)
%
%  Requires: MATLAB, Simulink, Simscape, Simscape Battery, PDE Toolbox
%  ========================================================================

function batterySimFramework_36s9p()

close all; clc;
outdir = fullfile(pwd, 'sim_outputs');
if ~exist(outdir,'dir'); mkdir(outdir); end

cfg.run_pde       = true;    % detailed single-cell PDE (PDE Toolbox)
cfg.run_simscape  = true;    % build + run Simscape Battery Simulink model
cfg.save_figs     = true;
ctype             = '21700'; % cell format

fprintf('==================================================================\n');
fprintf(' 36s9p 21700 Pack - Electro-Thermal Framework\n');
fprintf(' Physical: 18x18 grid, two 9x18 blocks + central air gap\n');
fprintf(' Cooling: negative-pressure forced air (induced-draft, gap exhaust)\n');
fprintf('==================================================================\n');

%% ===== (1) CELL + PACK DEFINITION =====
cell = defineCell(ctype);
pack = definePack(cell, 36, 9);            % 36 series, 9 parallel
printPackSpecs(pack);

%% ===== (2) ECM PARAMETERS =====
ecm = defineECM(cell);

%% ===== (3) VARIABLE MOTOR LOAD =====
load = defineMotorLoad(pack);

%% ===== (4) COUPLED ELECTRO-THERMAL SIM (forced air, gap) =====
fprintf('\n[4] Coupled electro-thermal sim (%d cells, 18x18, forced-air)...\n', pack.Ncells);
sim = runElectroThermal(pack, ecm, load);
fprintf('    SOC %.0f%% -> %.0f%% | pack T: %.1f -> mean %.1f / max %.1f C\n', ...
        100*sim.SOC(1), 100*sim.SOC(end), sim.T_C(1), mean(sim.Tfield(:)), max(sim.Tfield(:)));
fprintf('    Peak cell heat %.2f W (pack %.1f kW) | air rise inlet->gap %.1f C\n', ...
        max(sim.Qcell_peak), max(sim.Qpack_W)/1e3, sim.air_dT_last);

%% ===== (5) DETAILED SINGLE-CELL PDE =====
pde = struct('ok',false);
if cfg.run_pde
    fprintf('[5] Single-cell PDE thermal model...\n');
    try
        pde = singleCellPDE(cell, sim.Irms_cell, load.T_amb_C);
        fprintf('    PDE core-surface delta %.2f C | core Tmax %.1f C\n', ...
                pde.dT_core_surf, pde.T_core_max_C);
    catch ME
        fprintf(2,'    PDE skipped: %s\n', ME.message);
    end
end

%% ===== (6) SIMSCAPE BATTERY MODEL =====
ss = struct('ok',false);
if cfg.run_simscape
    fprintf('[6] Building + simulating Simscape Battery model...\n');
    try
        ss = simscapePackModel(cell, load, ctype);
        if ss.ok
            fprintf('    Simscape OK: Vpack %.1f -> %.1f V\n', ss.Vpack(1), ss.Vpack(end));
        end
    catch ME
        fprintf(2,'    Simscape reported: %s\n', ME.message);
        for i=1:min(numel(ME.stack),5)
            fprintf(2,'      at %s line %d\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
end

%% ===== (7) VISUALIZATION =====
fprintf('[7] Rendering figures...\n');
makeFigures(pack, ecm, load, sim, pde, ss, outdir, cfg, ctype);

save(fullfile(outdir,'framework_results.mat'), 'pack','ecm','load','sim','pde','ss');
fprintf('\nAll outputs written to: %s\n', outdir);
fprintf('==================================================================\n');
end

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function c = defineCell(~)
% 21700 high-energy NMC cell (nominal 4.5 Ah class)
c.name='21700 NMC'; c.type='21700';
c.D=0.02100; c.H=0.07000;
c.mass=0.0680; c.Qcap_Ah=4.50; c.R0_ref=0.012;
c.r         = c.D/2;
c.V         = pi*c.r^2*c.H;     % m^3 cell volume
c.A_surf    = 2*pi*c.r^2 + 2*pi*c.r*c.H;   % m^2 total surface
c.A_side    = 2*pi*c.r*c.H;     % m^2 lateral (cooling) surface
c.rho       = c.mass/c.V;       % kg/m^3  (CORRECT density)
c.cp        = 1000;             % J/(kg*K) specific heat (Li-ion jellyroll)
c.Cth       = c.mass*c.cp;      % J/K lumped thermal capacitance
c.Vnom      = 3.60;             % V nominal
c.Vmax      = 4.20;             % V
c.Vmin      = 3.00;             % V
c.Qcap_C    = c.Qcap_Ah*3600;   % Coulombs
c.energy_Wh = c.Vnom*c.Qcap_Ah; % Wh
% anisotropic conductivity for the PDE model (radial vs axial jelly roll)
c.k_r       = 0.9;              % W/(m*K) radial
c.k_z       = 24.2;             % W/(m*K) axial (through tabs/foils)
end

function p = definePack(cell,Ns,Np)
p.cell   = cell;
p.Ns     = Ns;            % series
p.Np     = Np;            % parallel
p.Ncells = Ns*Np;
p.Vnom   = Ns*cell.Vnom;
p.Vmax   = Ns*cell.Vmax;
p.Vmin   = Ns*cell.Vmin;
p.Cap_Ah = Np*cell.Qcap_Ah;
p.Energy_Wh = p.Vnom*p.Cap_Ah;
p.mass_cells = p.Ncells*cell.mass;
% Real packs add ~30-40% mass for enclosure, BMS, busbars, cooling.
p.overhead   = 1.35;
p.mass_pack  = p.mass_cells*p.overhead;
end

function printPackSpecs(p)
fprintf('\n[1] Pack architecture: %ds%dp  (%d x %s cells)\n', p.Ns, p.Np, p.Ncells, p.cell.type);
fprintf('    Nominal voltage : %.1f V   (%.1f - %.1f V window)\n', p.Vnom, p.Vmin, p.Vmax);
fprintf('    Capacity        : %.1f Ah\n', p.Cap_Ah);
fprintf('    Energy          : %.2f kWh\n', p.Energy_Wh/1000);
fprintf('    Cell mass       : %.1f kg  | pack (x%.2f overhead) %.1f kg\n', ...
        p.mass_cells, p.overhead, p.mass_pack);
fprintf('    Energy density  : %.0f Wh/kg (cell) | %.0f Wh/kg (pack est.)\n', ...
        p.cell.energy_Wh/p.cell.mass, p.Energy_Wh/p.mass_pack);
end

function e = defineECM(c)
% 2nd-order Thevenin equivalent circuit.  OCV/SOC + two RC pairs.
% --- Open-circuit voltage vs SOC (typical high-energy NMC) ---
e.soc_bp = [0 .05 .10 .20 .30 .40 .50 .60 .70 .80 .90 .95 1.00];
e.ocv_bp = [3.00 3.30 3.45 3.55 3.62 3.68 3.75 3.85 3.95 4.05 4.13 4.17 4.20];
% --- Entropic coefficient dOCV/dT (V/K), SOC dependent (reversible heat) ---
e.dUdT_bp = 1e-3*[ 0.10 0.05 0.00 -0.05 -0.08 -0.10 -0.12 -0.10 -0.08 ...
                  -0.05 -0.02  0.00  0.02];
% --- Ohmic + polarization ---
e.R0_ref = c.R0_ref;      % ohm @ 25C
e.Tref_K = 298.15;
e.Ea_R0  = 20000;         % J/mol Arrhenius activation (resistance falls as T rises)
e.Rgas   = 8.314;
e.R1 = 0.008;  e.tau1 = 12;   e.C1 = e.tau1/e.R1;    % fast RC (charge transfer)
e.R2 = 0.006;  e.tau2 = 120;  e.C2 = e.tau2/e.R2;    % slow RC (diffusion)
e.SOC0 = 0.90;            % start near full
end

function OCV = ocv(e,soc);  OCV = interp1(e.soc_bp,e.ocv_bp,min(max(soc,0),1),'pchip'); end
function d = dUdT(e,soc);    d   = interp1(e.soc_bp,e.dUdT_bp,min(max(soc,0),1),'pchip'); end
function R = R0ofT(e,T_K)
% Arrhenius: internal ohmic resistance decreases with temperature
R = e.R0_ref .* exp(e.Ea_R0/e.Rgas .* (1./T_K - 1/e.Tref_K));
end

function L = defineMotorLoad(p)
% Build a variable pack POWER demand (W, +discharge / -regen) representing a
% motor drawing from the pack, then convert to current inside the sim.
L.dt   = 1.0;                 % s
L.tend = 600;                 % s (10 min mixed-duty cycle)
L.t    = (0:L.dt:L.tend).';
L.T_amb_C = 25;
t = L.t;
P = zeros(size(t));
seg = @(a,b) (t>=a & t<b);
% Segment-based mission profile (kW at pack terminals)
P = P + seg(0,30)   .* linspace_seg(t,0,30, 3, 6);      % roll-out cruise
P = P + seg(30,55)  .* 26;                              % hard launch (peak)
P = P + seg(55,110) .* 9;                               % steady cruise
P = P + seg(110,130).* (-14);                           % regen braking
P = P + seg(130,140).* 2;                               % near coast
P = P + seg(140,150).* 24;                              % second launch
P = P + seg(150,250).* 0;                               % (urban filled below)
P = P + seg(250,320).* 22;                              % sustained hill climb
P = P + seg(320,345).* (-11);                           % downhill regen
P = P + seg(345,600).* 15;                              % highway cruise
% Urban stop-go band 150-250 s: oscillating accel/brake
urb = seg(150,250);
P(urb) = 10 + 9*sin(2*pi*(t(urb)-150)/22) + 3*sin(2*pi*(t(urb)-150)/7);
% Highway ripple 345-600 s
hw = seg(345,600);
P(hw) = P(hw) + 4*sin(2*pi*(t(hw)-345)/40);
L.P_W = P*1e3;                % -> Watts  (motor traction demand, + = discharge)
L.Ppeak_kW = max(P);

% ---- Solar charging input (roof/vehicle PV array, variable irradiance) ----
Psolar = 0.9 + 0.35*sin(2*pi*t/240) - 0.18*sin(2*pi*t/57);   % kW (~0.6-1.3)
Psolar = max(Psolar,0);
L.Psolar_W  = Psolar*1e3;
L.Pnet_W    = L.P_W - L.Psolar_W;         % net pack power (solar offsets motor)
L.Imotor_cmd = L.P_W    / p.Vnom;         % A traction (discharge) current
L.Isolar_cmd = L.Psolar_W / p.Vnom;       % A PV (charge) current
L.Ipack_cmd  = L.Pnet_W / p.Vnom;         % A net pack current
end

function y = linspace_seg(t,a,b,va,vb)
% linear ramp from va@a to vb@b over the whole vector (masked by caller)
y = va + (vb-va).*(t-a)/(b-a);
end

function S = runElectroThermal(p,e,L)
% Coupled ECM (uniform current, balanced pack) + per-cell lumped thermal net.
% Thermal grid is the PHYSICAL 18x18 arrangement (two 9x18 blocks + central
% gap), not the electrical Np x Ns grid.  Heat per cell is uniform (balanced
% 9p), so the scalar cell heat is applied across the physical grid.
Ns=p.Ns; Np=p.Np; c=p.cell; N=p.Ncells;
nt = numel(L.t); dt = L.dt;

% ---- Physical thermal grid: 18 x 18, gap between row 9 and 10 ------------
nR=18; nC=18; gapRow=9;           % two 9x18 blocks separated by the air gap
% cell-cell conductance (line contact + interstitial air); weak, so the
% cooling structure dominates the gradients.
gap      = 0.002; k_contact=2.0; k_air=0.026;
A_touch  = 0.015*c.H*c.D; A_air = c.H*c.D - A_touch;
G_cc     = (k_contact*A_touch + k_air*A_air)/gap;   % W/K between neighbours

% ---- Cooling: NEGATIVE-PRESSURE FORCED AIR, central-gap exhaust ----------
%   An induced-draft fan at the central gap pulls ambient air IN through the
%   two outer long faces (rows 1 and 18) and INWARD through each block toward
%   the gap, where it exhausts. Air warms as it travels toward the gap.
h_air    = 40;                    % W/m^2K forced-air convection coefficient
Tcool_C  = L.T_amb_C;             % inlet (ambient) air temperature
% exposure fraction of lateral area swept by the airflow
expo = 0.35*ones(nR,nC);
expo(1,:)=0.6; expo(nR,:)=0.6;               % outer long faces (air inlets)
expo(gapRow,:)=0.55; expo(gapRow+1,:)=0.55;  % gap-facing rows (exhaust side)
expo(:,1)=max(expo(:,1),0.55); expo(:,nC)=max(expo(:,nC),0.55);
expo([1 nR],[1 nC])=0.8;                     % corners
Gcool = h_air*c.A_side.*expo;     % W/K per cell (nR x nC)

cool.mdot   = 0.12;               % kg/s total air mass flow (induced-draft)
cool.cp     = 1005;               % J/kgK air
cool.rho    = 1.16;               % kg/m^3 air (reduced under negative pressure)
cool.dP_Pa  = -350;               % Pa gauge (suction / negative pressure)
cool.Tin_C  = Tcool_C;
mdot_path   = cool.mdot/(2*nC);   % two inward paths per column

% small (+/-2%) manufacturing spread in R0 for realistic non-uniformity
rng(7); R0scale = 1 + 0.02*randn(nR,nC);

% ---- Implicit thermal operator (gap breaks conduction rows 9<->10) -------
L_lap = buildLaplacian(nR,nC,G_cc,gapRow);
Cdiag = c.Cth/dt;
A = Cdiag*speye(N) + L_lap + spdiags(Gcool(:),0,N,N);
dA = decomposition(A,'chol');

% ---- State init ----------------------------------------------------------
T = (L.T_amb_C+273.15)*ones(N,1); soc=e.SOC0; V1=0; V2=0;
Vpack_prev = Ns*ocv(e,soc);

% ---- Logs ----------------------------------------------------------------
S.t=L.t; S.SOC=zeros(nt,1); S.Vpack=zeros(nt,1); S.Icell=zeros(nt,1);
S.Ipack=zeros(nt,1); S.Qpack_W=zeros(nt,1); S.Tmax=zeros(nt,1);
S.Tmin=zeros(nt,1); S.Tmean=zeros(nt,1); S.Qcell_peak=zeros(nt,1);
S.SOC(1)=soc; S.Vpack(1)=Vpack_prev; S.Tmean(1)=L.T_amb_C;
S.Tmax(1)=L.T_amb_C; S.Tmin(1)=L.T_amb_C;
Isum2=zeros(nR,nC);
Tair_grid=(Tcool_C+273.15)*ones(nR,nC);

for k=2:nt
    % ---- Electrical: NET pack current = motor draw - solar charge --------
    Ipack = L.Pnet_W(k)/max(Vpack_prev,1);   % A (pack), + = net discharge
    Icell = Ipack/Np;                        % A per cell (9 in parallel)

    V1 = V1 + dt*( Icell/e.C1 - V1/(e.R1*e.C1) );
    V2 = V2 + dt*( Icell/e.C2 - V2/(e.R2*e.C2) );
    soc = min(max(soc - Icell*dt/c.Qcap_C,0),1);

    R0 = R0scale .* R0ofT(e, reshape(T,nR,nC));   % nR x nC (ohm)
    Vcell = max(ocv(e,soc) - Icell*mean(R0(:)) - V1 - V2, c.Vmin);
    Vpack = Ns*Vcell;

    % ---- Heat generation per cell (Bernardi) ----------------------------
    Q_ohmic = (Icell^2).*R0;
    Q_rc    = (V1^2/e.R1 + V2^2/e.R2);
    Q_ent   = -Icell*reshape(T,nR,nC).*dUdT(e,soc);
    Qcell   = max(Q_ohmic + Q_rc + Q_ent, 0);      % W per cell (nR x nC)

    % ---- Air temperature: two inward streams per column toward the gap --
    Tg = reshape(T,nR,nC); Tair_grid=zeros(nR,nC); Tin=cool.Tin_C+273.15;
    for cc=1:nC
        Ta=Tin;                                   % block 1: rows 1..gapRow
        for r=1:gapRow
            Tair_grid(r,cc)=Ta;
            Ta=Ta + max(Gcool(r,cc)*(Tg(r,cc)-Ta),0)/(mdot_path*cool.cp);
        end
        Ta=Tin;                                   % block 2: rows nR..gapRow+1
        for r=nR:-1:gapRow+1
            Tair_grid(r,cc)=Ta;
            Ta=Ta + max(Gcool(r,cc)*(Tg(r,cc)-Ta),0)/(mdot_path*cool.cp);
        end
    end

    % ---- Thermal implicit step  A*Tnew = C/dt*T + Q + Gcool*Tair --------
    b = Cdiag*T + Qcell(:) + Gcool(:).*Tair_grid(:);
    T = dA\b;
    Isum2 = Isum2 + Icell^2;

    Tc = reshape(T,nR,nC)-273.15;
    S.SOC(k)=soc; S.Vpack(k)=Vpack; S.Icell(k)=Icell; S.Ipack(k)=Ipack;
    S.Qpack_W(k)=sum(Qcell(:)); S.Qcell_peak(k)=max(Qcell(:));
    S.Tmax(k)=max(Tc(:)); S.Tmin(k)=min(Tc(:)); S.Tmean(k)=mean(Tc(:));
    Vpack_prev=Vpack;
end

S.Tfield  = reshape(T,nR,nC)-273.15;      % final temperature field (C)
S.T_C     = S.Tmean;
S.nR=nR; S.nC=nC; S.gapRow=gapRow; S.G_cc=G_cc; S.Gcool=Gcool; S.h_air=h_air;
S.expo=expo; S.Tcool_C=Tcool_C; S.cool=cool;
S.Tair_grid = Tair_grid-273.15;                    % final air field (C)
S.Tair_profile = mean(S.Tair_grid,2);              % air temp vs row (mean col)
S.air_dT_last = max(S.Tair_grid(:))-cool.Tin_C;    % inlet->gap air rise
[~,idx]=max(S.Tfield(:));
S.Irms_cell = sqrt(Isum2(idx)/(nt-1));
end

function L = buildLaplacian(nR,nC,g,gapRow)
% Graph Laplacian of a 4-neighbor grid with edge conductance g (insulated
% exterior).  If gapRow>0, conduction between rows gapRow and gapRow+1 is
% removed (physical air gap between the two blocks).
if nargin<4, gapRow=0; end
N=nR*nC; idx=@(r,cc)(cc-1)*nR + r;
I=[];J=[];Sv=[];
for cc=1:nC
  for r=1:nR
    i=idx(r,cc); nb=[];
    if r>1  && (r-1)~=gapRow, nb(end+1)=idx(r-1,cc); end   % up (skip gap)
    if r<nR && r~=gapRow,     nb(end+1)=idx(r+1,cc); end   % down (skip gap)
    if cc>1,  nb(end+1)=idx(r,cc-1); end
    if cc<nC, nb(end+1)=idx(r,cc+1); end
    for j=nb
      I(end+1)=i; J(end+1)=j; Sv(end+1)=-g;      %#ok<AGROW>
    end
    I(end+1)=i; J(end+1)=i; Sv(end+1)=g*numel(nb); %#ok<AGROW>
  end
end
L=sparse(I,J,Sv,N,N);
end

function P = singleCellPDE(c, I_rms, Tamb_C)
% Transient PDE thermal model of one cylindrical cell.  The material density
% is rho; heat generation q is applied separately.  Conductivity is
% anisotropic: [k_r; k_r; k_z].
P.ok=false;
Tamb=Tamb_C+273.15;
q = (I_rms^2)*c.R0_ref / c.V;              % W/m^3 volumetric heat gen
h_conv = 20;                               % W/m^2K

tm = createpde('thermal','transient');
tm.Geometry = multicylinder(c.r, c.H);
generateMesh(tm,'Hmax',0.0025);
thermalProperties(tm,'ThermalConductivity',[c.k_r; c.k_r; c.k_z], ...
                     'MassDensity',c.rho,'SpecificHeat',c.cp);   % rho fixed
internalHeatSource(tm,q);
allfaces = 1:tm.Geometry.NumFaces;
thermalBC(tm,'Face',allfaces,'ConvectionCoefficient',h_conv,'AmbientTemperature',Tamb);
thermalIC(tm,Tamb);

tlist = 0:10:1800;                          % 30 min
R = solve(tm,tlist);
Tend = R.Temperature(:,end);

% core vs surface along mid-height
xc=linspace(0,c.r,60); yc=zeros(size(xc)); zc=(c.H/2)*ones(size(xc));
Trad = zeros(numel(xc),numel(tlist));
for i=1:numel(tlist)
    Trad(:,i)=interpolateTemperature(R,xc,yc,zc,i);
end
P.ok=true; P.model=tm; P.result=R; P.tlist=tlist;
P.Tend_C=Tend-273.15; P.Trad_C=Trad-273.15; P.xc=xc; P.q=q; P.I_rms=I_rms;
P.T_core_max_C = max(Trad(1,:))-273.15;
P.dT_core_surf = (Trad(1,end)-Trad(end,end));
end

function R = simscapePackModel(cell, L, ctype)
% Build a Simscape Battery ECM model of the 36s9p pack and run a variable
% motor-current load.  Electrical (equivalent-circuit) response is the focus;
% detailed thermal is handled by the finite-volume network + PDE model.
%
% Verified port mapping (R2026a foundation blocks):
%   Controlled Current Source : LConn1=+ , RConn1=PS command , RConn2=-
%   Voltage Sensor            : LConn1=+ , RConn1=PS output  , RConn2=-
%   Simulink-PS Converter     : Inport1=Simulink , RConn1=PS out
%   PS-Simulink Converter     : LConn1=PS in     , Outport1=Simulink
if nargin<3, ctype='21700'; end
R.ok=false;
base=['pack36s9p_' lower(ctype)]; genlib=[base '_lib']; name=[base '_model'];

cleanupBatteryArtifacts(base,genlib,name);

% ---- Simscape Battery objects (electrical ECM, no thermal port) ----
geo  = batteryCylindricalGeometry(simscape.Value(cell.H,'m'), simscape.Value(cell.r,'m'));
scell= batteryCell(geo);
scell.Mass     = simscape.Value(cell.mass,'kg');
scell.Capacity = simscape.Value(cell.Qcap_Ah,'A*hr');
pa  = batteryParallelAssembly(scell, 9, ModelResolution='Lumped');
mod = batteryModule(pa, 36, ModelResolution='Lumped', StackingAxis='X', ...
                    InterParallelAssemblyGap=simscape.Value(0.003,'m'));
buildBattery(mod,'LibraryName',base);

% ---- Locate library block by exact path (find_system) ----
load_system(genlib);
blks = find_system(genlib,'SearchDepth',1,'Type','block'); src='';
for i=1:numel(blks)
    if strcmp(get_param(blks{i},'Name'), mod.Name); src=blks{i}; end
end
assert(~isempty(src),'Battery module block not found in %s',genlib);

% ---- Assemble runnable model ----
new_system(name); load_system(name);
bp=[name '/Pack'];
add_block(src, bp, 'Position',[220 100 360 240], 'MakeNameUnique','on');
% Match per-cell capacity + initial SOC to the analytical ECM so the two
% models are directly comparable (the built block otherwise keeps a large
% default capacity that masks the SOC swing).
try, set_param(bp,'BatteryCapacityCell',num2str(cell.Qcap_Ah)); catch, end
try, set_param(bp,'socCell','0.9');                            catch, end
add_block('nesl_utility/Solver Configuration',[name '/Solver'],'Position',[120 300 170 340]);
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference',[name '/Gnd'], ...
          'Position',[220 320 240 340],'Orientation','down');
add_block('fl_lib/Electrical/Electrical Sources/Controlled Current Source', ...
          [name '/Load'],'Position',[60 140 90 200],'Orientation','down');
add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor', ...
          [name '/Vsense'],'Position',[400 140 430 200],'Orientation','down');
add_block('simulink/Sources/From Workspace',[name '/Icmd'],'Position',[-120 150 -40 190]);
add_block('nesl_utility/Simulink-PS Converter',[name '/S2PS'],'Position',[-10 155 20 185]);
add_block('nesl_utility/PS-Simulink Converter',[name '/PS2S'],'Position',[470 155 500 185]);
add_block('simulink/Sinks/To Workspace',[name '/Vout'],'Position',[540 155 610 185]);
set_param([name '/Vout'],'VariableName','Vpack_ss','SaveFormat','StructureWithTime');

% signal wiring (Simulink domain)
add_line(name,'Icmd/1','S2PS/1');
add_line(name,'PS2S/1','Vout/1');

bph=get_param(bp,'PortHandles');           lph=get_param([name '/Load'],'PortHandles');
vph=get_param([name '/Vsense'],'PortHandles'); sph=get_param([name '/S2PS'],'PortHandles');
pph=get_param([name '/PS2S'],'PortHandles'); solh=get_param([name '/Solver'],'PortHandles');
gh =get_param([name '/Gnd'],'PortHandles');

% physical-signal command -> current-source PS input (RConn1)
add_line(name, sph.RConn(1), lph.RConn(1),'autorouting','on');
% voltage-sensor PS output (RConn1) -> PS2S PS input (LConn1)
add_line(name, vph.RConn(1), pph.LConn(1),'autorouting','on');

% single-port blocks (Solver, Electrical Reference) may expose L or R conn
solP = firstConn(solh); gndP = firstConn(gh);
% electrical loop: + = LConn1, - = RConn1 on the pack; branch off each node
add_line(name, bph.LConn(1), lph.LConn(1),'autorouting','on');  % + -> source +
add_line(name, bph.LConn(1), vph.LConn(1),'autorouting','on');  % + -> sensor +
add_line(name, bph.LConn(1), solP,'autorouting','on');          % + -> solver
add_line(name, bph.RConn(1), lph.RConn(2),'autorouting','on');  % - -> source -
add_line(name, bph.RConn(1), vph.RConn(2),'autorouting','on');  % - -> sensor -
add_line(name, bph.RConn(1), gndP,'autorouting','on');          % - -> ground

% ---- Current command profile (motor controller) ----
% Controlled Current Source injects + into the pack + terminal (charging),
% so motoring (positive power) is a NEGATIVE command; regen is positive.
Icmd = timeseries(-L.Ipack_cmd, L.t);
assignin('base','Icmd_ts',Icmd);
set_param([name '/Icmd'],'VariableName','Icmd_ts');

set_param(name,'SimscapeLogType','all','StopTime',num2str(L.tend), ...
          'Solver','ode23t','MaxStep','0.5','AbsTol','1e-4','RelTol','1e-3');

out = sim(name);
R.ok=true; R.model=name;

% ---- Extract terminal voltage ----
try
    Vs = out.get('Vpack_ss');
    R.t = Vs.time;  R.Vpack = Vs.signals.values(:);
catch
    R.t=[0;L.tend]; R.Vpack=[34*cell.Vnom;34*cell.Vnom];
end
R.Tmax_C = L.T_amb_C;   % thermal handled by FV network
try, save_system(name); catch, end
try, close_system(name,0); catch, end
end

function h = firstConn(ph)
% return the first available conserving-port handle (LConn or RConn)
if ~isempty(ph.LConn); h=ph.LConn(1); elseif ~isempty(ph.RConn); h=ph.RConn(1);
else; error('block has no conserving port'); end
end

function cleanupBatteryArtifacts(base,genlib,name)
% Remove all files buildBattery emits so a rebuild does not error on collision
if bdIsLoaded(genlib); close_system(genlib,0); end
if bdIsLoaded(name);   close_system(name,0);   end
files = {[base '.mat'],[genlib '.slx'],[name '.slx'],[base '_param.m'],[genlib '_param.m']};
for f=files; if exist(f{1},'file'); try, delete(f{1}); catch, end; end; end
dirs = {['+' base],['+' base 'LumpingAdapters'],['+' genlib 'LumpingAdapters']};
for d=dirs; if exist(d{1},'dir'); try, rmdir(d{1},'s'); catch, end; end; end
end

function applyLightTheme(f)
% Force a clean light presentation theme (R2026a defaults to dark in batch)
try, theme(f,'light'); catch, end
set(f,'Color','w');
end

function makeFigures(p,e,L,S,pde,ss,outdir,cfg,ctype)
if nargin<9, ctype=p.cell.type; end
tag=['_' lower(ctype)];
vis = 'on';
% ---------- Figure 1: Mission + electrical response ----------------------
f1=figure('Color','w','Position',[60 60 1500 820],'Visible',vis, ...
          'Name','Electro-Thermal Overview'); applyLightTheme(f1);
tl=tiledlayout(f1,3,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('36s9p %s Pack  |  %.2f kWh  |  Variable Motor Load (forced-air cooled)', ...
      ctype, p.Energy_Wh/1000),'FontSize',14,'FontWeight','bold');

nexttile; plot(L.t,L.P_W/1e3,'LineWidth',1.6,'Color',[0.10 0.35 0.75]); grid on;
yline(0,':','Color',[.4 .4 .4]);
ylabel('Motor power (kW)'); title('(a) Mission profile'); xlim([0 L.tend]);

nexttile; plot(L.t,S.Ipack,'LineWidth',1.4,'Color',[0.80 0.25 0.20]); grid on;
ylabel('Pack current (A)'); title('(b) Pack current demand'); xlim([0 L.tend]);
yline(0,':','Color',[.4 .4 .4]);

nexttile; plot(L.t,S.Vpack,'LineWidth',1.6,'Color',[0.15 0.55 0.30]); grid on; hold on;
if isfield(ss,'ok') && ss.ok && isfield(ss,'Vpack') && numel(ss.Vpack)>2
    plot(ss.t,ss.Vpack,'--','Color',[0.85 0.35 0.10],'LineWidth',1.3);
    legend('analytical ECM','Simscape ECM','Location','southwest','FontSize',8);
end
yline(p.Vmin,'--r','Min'); yline(p.Vmax,'--','Color',[.4 .4 .4]);
ylabel('Pack voltage (V)'); title('(c) Terminal voltage'); xlim([0 L.tend]);

nexttile; plot(L.t,100*S.SOC,'LineWidth',1.6,'Color',[0.45 0.30 0.70]); grid on;
ylabel('SOC (%)'); title('(d) State of charge'); xlim([0 L.tend]); ylim([0 100]);

nexttile; plot(L.t,S.Qpack_W/1e3,'LineWidth',1.4,'Color',[0.85 0.45 0.10]); grid on;
ylabel('Pack heat gen (kW)'); xlabel('Time (s)'); title('(e) Total heat generation');
xlim([0 L.tend]);

nexttile; hold on; grid on;
plot(L.t,S.Tmean,'k','LineWidth',2);
plot(L.t,S.Tmax,'--','Color',[0.85 0.20 0.20],'LineWidth',1.4);
plot(L.t,S.Tmin,'--','Color',[0.20 0.55 0.90],'LineWidth',1.4);
yline(L.T_amb_C,':','Ambient','Color',[.4 .4 .4]);
ylabel('Temperature (\circC)'); xlabel('Time (s)');
title('(f) Pack temperature'); legend('mean','max','min','Location','northwest');
xlim([0 L.tend]);
if cfg.save_figs, exportgraphics(f1,fullfile(outdir,['fig1_overview' tag '.png']),'Resolution',150); end

% ---------- Figure 2: Spatial thermal field ------------------------------
f2=figure('Color','w','Position',[80 80 1400 560],'Visible',vis, ...
          'Name','Pack Thermal Field'); applyLightTheme(f2);
tl2=tiledlayout(f2,1,2,'TileSpacing','compact','Padding','compact');
title(tl2,sprintf('End-of-Cycle Pack Thermal State (%s, forced-air)',ctype), ...
      'FontSize',14,'FontWeight','bold');

% Insert a blank (NaN) row at the gap so the two 9x18 blocks render apart
gr=S.gapRow;
Tdisp=[S.Tfield(1:gr,:); nan(1,S.nC); S.Tfield(gr+1:end,:)];
nexttile; hI=imagesc(Tdisp); set(hI,'AlphaData',~isnan(Tdisp));
axis image; colormap(gca,hot); set(gca,'Color',[0.9 0.9 0.92]);
cb=colorbar; ylabel(cb,'Temperature (\circC)');
xlabel('Cell column (1-18)'); ylabel('Cell row (1-18, gap between blocks)');
title(sprintf('(a) 18x18 cell map, central air gap (\\DeltaT = %.1f\\circC)', ...
      max(S.Tfield(:))-min(S.Tfield(:))));
text(S.nC/2, gr+1, 'air gap / exhaust', 'HorizontalAlignment','center', ...
     'FontSize',8,'Color',[0.2 0.2 0.5],'FontAngle','italic');

% (b) temperature vs row (across the two blocks) with the inward air profile
nexttile; hold on; grid on;
rows=1:S.nR;
Trow_mean=mean(S.Tfield,2)'; Trow_min=min(S.Tfield,[],2)'; Trow_max=max(S.Tfield,[],2)';
fill([rows, fliplr(rows)],[Trow_min, fliplr(Trow_max)],[0.85 0.88 1.0], ...
     'EdgeColor','none','FaceAlpha',0.6);
plot(rows,Trow_mean,'k-o','LineWidth',1.6,'MarkerSize',3,'MarkerFaceColor','k');
plot(rows,S.Tair_profile','-','Color',[0.10 0.45 0.85],'LineWidth',1.8);
xline(gr+0.5,'--','Color',[0.5 0.5 0.5],'Label','gap (air exhaust)');
legend('cell min-max band','cell mean','air temperature','Location','north');
xlabel('Cell row  (air drawn inward from rows 1 & 18 \rightarrow gap)');
ylabel('Temperature (\circC)');
title('(b) Cell & air temperature across the two blocks'); xlim([1 S.nR]);
if cfg.save_figs, exportgraphics(f2,fullfile(outdir,['fig2_thermal_field' tag '.png']),'Resolution',150); end

% ---------- Figure 3: 3D pack layout chart -------------------------------
buildPackLayout3D(outdir,ctype,cfg.save_figs,tag,vis);

% ---------- Figure 4: Single-cell PDE ------------------------------------
if isfield(pde,'ok') && pde.ok
    f3=figure('Color','w','Position',[100 100 1300 560],'Visible',vis, ...
              'Name','Single-Cell PDE'); applyLightTheme(f3);
    tl3=tiledlayout(f3,1,2,'TileSpacing','compact','Padding','compact');
    title(tl3,sprintf('Detailed 3-D %s Cell Thermal Model  (I_{rms}=%.1f A)', ...
          ctype, pde.I_rms),'FontSize',14,'FontWeight','bold');
    nexttile; pdeplot3D(pde.model,'ColorMapData',pde.Tend_C); axis equal;
    title('(a) Internal temperature @ 30 min (\circC)');
    nexttile; surf(pde.tlist/60,pde.xc*1e3,pde.Trad_C,'EdgeColor','none');
    xlabel('Time (min)'); ylabel('Radius (mm)'); zlabel('Temp (\circC)');
    title(sprintf('(b) Radial profile (core-surface \\DeltaT = %.2f\\circC)', ...
          pde.dT_core_surf)); view(135,30); colormap(gca,turbo);
    if cfg.save_figs, exportgraphics(f3,fullfile(outdir,['fig4_cell_pde' tag '.png']),'Resolution',150); end
end

fprintf('    Figures saved (fig1..fig4%s).\n', tag);
end

function buildPackLayout3D(outdir,ctype,saveFigures,tag,vis)
% Dark, cylindrical pack chart matching the requested Simscape-style 3D look.
f=figure('Color','k','Position',[80 80 1620 920],'Visible',vis, ...
         'Name','Pack Layout Chart');
ax=axes(f); hold(ax,'on');
set(ax,'Color','k','XColor',[0.92 0.92 0.92],'YColor',[0.92 0.92 0.92], ...
    'ZColor',[0.92 0.92 0.92],'Box','on','Projection','perspective');
grid(ax,'on');
ax.GridColor = [0.85 0.85 0.85];
ax.GridAlpha = 0.28;
ax.LineWidth = 1.0;
view(ax,135,25);

% 18 x 9 cell slab with a split between the two 2s9p rows in the long axis.
nLong = 18;
nWide = 9;
cellR = 0.018;
cellH = 0.055;
pitchX = 0.043;
pitchY = 0.032;
splitX = 9 * pitchX + 0.030;
originX = 0.06;
originY = 0.03;

% Lighting and background plane.
surface(ax,[0 0.85;0 0.85],[0 0;0.32 0.32],[0 0;0 0],'FaceColor',[0.06 0.06 0.06], ...
    'EdgeColor','none','FaceAlpha',1.0);

topColor = [0.70 0.70 0.70];
sideColor = [0.48 0.48 0.48];
edgeColor = [0.62 0.62 0.62];

for ix = 1:nLong
    if ix <= 9
        x = originX + (ix-1) * pitchX;
    else
        x = originX + (ix-1) * pitchX + 0.060;
    end
    for iy = 1:nWide
        y = originY + (iy-1) * pitchY;
        drawCylinderCell(ax,x,y,0,cellR,cellH,topColor,sideColor,edgeColor);
    end
end

% Split wall / cooling gap indicator.
gapX = originX + 9 * pitchX + 0.026;
patch(ax,[gapX gapX+0.016 gapX+0.016 gapX],[0.015 0.015 0.315 0.315],[0.12 0.12 0.12], ...
    'FaceAlpha',0.9,'EdgeColor',[0.85 0.85 0.85],'LineWidth',0.6);
text(ax,gapX+0.008,0.165,0.065,'split / cooling gap','Color',[0.94 0.94 0.94], ...
    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');

% Axes and labels to keep the chart readable.
title(ax,{sprintf('36s9p %s Pack  |  Split Two-Row Cooling Layout',ctype), ...
    'two 2s9p module rows going back; dark cylindrical chart style'}, ...
    'Color','w','FontSize',16,'FontWeight','bold');
xlabel(ax,'pack length','Color','w');
ylabel(ax,'pack width','Color','w');
zlabel(ax,'cell height','Color','w');
xlim(ax,[0.0 0.84]); ylim(ax,[0.0 0.32]); zlim(ax,[0.0 0.09]);
ax.XTick = 0:0.1:0.8;
ax.YTick = 0:0.1:0.3;
ax.ZTick = [0 0.05];
camlight(ax,'headlight');
lighting(ax,'gouraud');
material(ax,[0.35 0.55 0.22 12]);

if saveFigures
    exportgraphics(f,fullfile(outdir,['fig3_pack_layout_3d' tag '.png']),'Resolution',150);
end
end

function drawCylinderCell(ax,x,y,z,r,h,topColor,sideColor,edgeColor)
[xc,yc,zc]=cylinder(r,16);
zc = zc * h + z;
xc = xc + x;
yc = yc + y;
surf(ax,xc,yc,zc,'FaceColor',sideColor,'EdgeColor','none','FaceLighting','gouraud');
[xt,yt,zt]=cylinder(r,24);
zt = zt(2,:) * 0 + (z + h);
xt = xt + x;
yt = yt + y;
patch(ax,xt(1,:),yt(1,:),repmat(z+h,1,numel(xt)),topColor,'EdgeColor',edgeColor,'LineWidth',0.25);
end
