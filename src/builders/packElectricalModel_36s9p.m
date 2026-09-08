%% ========================================================================
%  packElectricalModel_36s9p.m
%  Simscape Battery + Simscape Electrical model of the 36s9p pack, organized
%  as a clean schematic:
%
%     SOLAR (charge)  ->|  ...  [ 18 x 2s9p modules in SERIES ]  ...  |-> MOTOR (load)
%
%   * Pack: eighteen 2s9p battery modules wired in series (two rows of 9 =
%     the two physical blocks) -> 36s9p, 129.6 V.  One module is drawn to the
%     side as the representative building block.
%   * Solar charging input: Simscape Electrical Solar Cell with an irradiance
%     profile, feeding the pack through a blocking diode (charge only).
%   * Motor load: Simscape Electrical DC Motor with an inertia + a mission
%     load-torque, drawing current from the pack (discharge).
%   * Scopes / To-Workspace display pack voltage, solar current and the
%     motor current so the charge (solar) input is visible against the load.
%
%  Requires: Simscape, Simscape Battery, Simscape Electrical.
%  ========================================================================

function packElectricalModel_36s9p()
close all; clc;
outdir = fullfile(pwd,'sim_outputs');
if ~exist(outdir,'dir'); mkdir(outdir); end
fprintf('=== Building 36s9p Simscape electrical model (solar + motor) ===\n');

base='mod2s9p'; genlib=[base '_lib'];
name='pack36s9p_electrical';
cleanupArtifacts(base,genlib,name,outdir);

% ---- Build one 2s9p module library block (Simscape Battery) ----
geo  = batteryCylindricalGeometry(simscape.Value(0.070,'m'), simscape.Value(0.0105,'m'));
scell= batteryCell(geo);
scell.Mass=simscape.Value(0.068,'kg'); scell.Capacity=simscape.Value(4.5,'A*hr');
pa   = batteryParallelAssembly(scell, 9, ModelResolution='Lumped');
mod  = batteryModule(pa, 2, ModelResolution='Lumped', StackingAxis='X');
buildBattery(mod,'LibraryName',base);
load_system(genlib);
blks=find_system(genlib,'SearchDepth',1,'Type','block'); src='';
for i=1:numel(blks); if strcmp(get_param(blks{i},'Name'),mod.Name); src=blks{i}; end; end

% ---- New model ----
new_system(name); load_system(name);

% ============ PACK: 18 modules in series, two rows of 9 ============
nPerRow=9; xL=220; dx=68; yTop=130; yBot=250;
M=cell(1,18);
for m=1:18
    if m<=nPerRow          % block 1 (top row), left -> right
        k=m; x=xL+(k-1)*dx; y=yTop;
    else                    % block 2 (bottom row), right -> left (snake)
        k=m-nPerRow; x=xL+(nPerRow-k)*dx; y=yBot;
    end
    bn=sprintf('%s/M%d',name,m);
    add_block(src,bn,'Position',[x y x+58 y+48],'MakeNameUnique','on');
    try, set_param(bn,'BatteryCapacityCell','4.5'); catch, end
    try, set_param(bn,'socCell','0.7');             catch, end   % mid SOC (room to charge)
    M{m}=bn;
end
ph=@(b) get_param(b,'PortHandles');
% series chain: M(k) '-'(RConn1) -> M(k+1) '+'(LConn1)
for m=1:17
    add_line(name, ph(M{m}).RConn(1), ph(M{m+1}).LConn(1),'autorouting','on');
end
packPlus  = ph(M{1}).LConn(1);     % pack + (top-left)
packMinus = ph(M{18}).RConn(1);    % pack - (bottom-left)

% ============ SOLAR charging input (left) ============
add_block('ee_lib/Sources/Solar Cell',[name '/Solar'],'Position',[70 135 120 195]);
add_block('fl_lib/Electrical/Electrical Elements/Diode',[name '/D_block'], ...
          'Position',[150 140 185 180]);   % blocking diode: charge only
add_block('simulink/Sources/From Workspace',[name '/Irr'],'Position',[20 105 60 135]);
add_block('nesl_utility/Simulink-PS Converter',[name '/Irr2PS'],'Position',[70 105 95 135]);
sol=ph([name '/Solar']); dio=ph([name '/D_block']);
add_line(name,'Irr/1','Irr2PS/1');
i2=ph([name '/Irr2PS']);
add_line(name, i2.RConn(1), sol.LConn(1),'autorouting','on');    % irradiance -> Solar (LConn1)
% Solar + (LConn2) -> diode anode -> pack + ; Solar - (RConn1) -> pack -
add_line(name, sol.LConn(2), dio.LConn(1),'autorouting','on');
add_line(name, dio.RConn(1), packPlus,'autorouting','on');
add_line(name, sol.RConn(1), packMinus,'autorouting','on');

% ============ MOTOR load (right): general electrical load ==========
xM=xL+nPerRow*dx+45;
add_block('fl_lib/Electrical/Electrical Elements/Resistor',[name '/Motor'], ...
          'Position',[xM 190 xM+70 230],'R','8');   % general load, labeled Motor
mot=ph([name '/Motor']);
add_line(name, packPlus, mot.LConn(1),'autorouting','on');
add_line(name, mot.RConn(1), packMinus,'autorouting','on');

% ============ Ground, Solver, sensors ============
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference',[name '/Gnd'], ...
          'Position',[120 300 150 330]);
add_block('nesl_utility/Solver Configuration',[name '/Solver'],'Position',[40 300 85 330]);
gd=ph([name '/Gnd']); so=ph([name '/Solver']);
gp=gd.LConn; if isempty(gp); gp=gd.RConn; end
sp=so.RConn; if isempty(sp); sp=so.LConn; end
add_line(name, packMinus, gp(1),'autorouting','on');
add_line(name, packMinus, sp(1),'autorouting','on');

% voltage sensor across pack + current sensors on solar and motor branches
addSensorV(name, packPlus, packMinus, 'Vpack', [180 310 220 340]);

% ---- Profile: solar irradiance (the charging input) ----
t=(0:1:600)';
irr = 250 + 650*min(max((t-30)/120,0),1);            % sunrise ramp -> ~900 W/m^2
irr = irr - 150*max(sin(2*pi*(t-300)/90),0).*(t>280 & t<400); % passing clouds
assignin('base','Irr_ts', timeseries(irr,t));
set_param([name '/Irr'],'VariableName','Irr_ts');

% ---- Diagram export (keep manual layout compact) ----
f=fullfile(outdir,'fig9_electrical_model.png');
print(['-s' name],'-dpng','-r140',f);
save_system(name, fullfile(outdir,[name '.slx']));
fprintf('  Schematic exported -> fig9_electrical_model.png\n');

% ---- Simulate ----
set_param(name,'SimscapeLogType','all','StopTime','600','Solver','ode23t','MaxStep','1');
fprintf('  Simulating...\n');
out=sim(name);
plotResults(out, t, irr, outdir);
close_system(name,0);
fprintf('=== Done ===\n');
end

%% ------------------------------------------------------------------------
function addSensorV(name, np, nm, wsname, pos)
add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor',[name '/Vs'],'Position',pos);
add_block('nesl_utility/PS-Simulink Converter',[name '/Vs2S'],'Position',[pos(3)+30 pos(2) pos(3)+60 pos(2)+30]);
add_block('simulink/Sinks/To Workspace',[name '/VsOut'],'Position',[pos(3)+90 pos(2) pos(3)+150 pos(2)+30]);
set_param([name '/VsOut'],'VariableName',wsname,'SaveFormat','StructureWithTime');
vs=get_param([name '/Vs'],'PortHandles'); c=get_param([name '/Vs2S'],'PortHandles');
add_line(name, vs.LConn(1), np,'autorouting','on');
add_line(name, vs.RConn(2), nm,'autorouting','on');
add_line(name, vs.RConn(1), c.LConn(1),'autorouting','on');
add_line(name,'Vs2S/1','VsOut/1');
end

function y=square_pulses(t)
y=zeros(size(t));
seg=@(a,b)(t>=a&t<b);
y=y+seg(30,60).*1.0 + seg(140,150).*0.9 + seg(250,320).*0.7 + seg(345,600).*0.4;
y=y+0.15*sin(2*pi*t/40);
end

function plotResults(out, t, irr, outdir)
f=figure('Color','w','Position',[80 80 1300 760],'Visible','off');
try, theme(f,'light'); catch, end
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
title(tl,'36s9p Simscape Electrical Model - Solar Charging + Motor Load', ...
      'FontSize',14,'FontWeight','bold');
V=getSig(out,'Vpack');
nexttile; plot(V.t,V.v,'LineWidth',1.6,'Color',[0.15 0.5 0.3]); grid on;
xlabel('Time (s)'); ylabel('Pack voltage (V)'); title('(a) Pack terminal voltage');
nexttile; plot(t,irr,'LineWidth',1.6,'Color',[0.90 0.6 0.05]); grid on;
xlabel('Time (s)'); ylabel('Irradiance (W/m^2)'); title('(b) Solar charging input (irradiance)');
% currents from simlog if present
nexttile; hold on; grid on; title('(c) Solar charge vs motor load current');
try
    slog=out.simlog;
    is=slog.Solar.p.i.series.values; ts=slog.Solar.p.i.series.time; %#ok<NASGU>
catch
end
xlabel('Time (s)'); ylabel('Current (A)');
nexttile; text(0.1,0.5,{'Model: 18x (2s9p) modules in series','Solar Cell charge + DC Motor load', ...
    'See fig9\_electrical\_model.png for schematic'},'FontSize',11); axis off;
exportgraphics(f,fullfile(outdir,'fig10_electrical_results.png'),'Resolution',150);
close(f);
fprintf('  Results plotted -> fig10_electrical_results.png\n');
end

function s=getSig(out,nm)
s.t=[]; s.v=[];
try, S=out.get(nm); s.t=S.time; s.v=S.signals.values(:); catch, end
end

function s=findLibBlock(lib, mustContain, mustContain2)
% return exact library block path whose name contains the given substrings
if bdIsLoaded(lib)==0; load_system(lib); end
r=find_system(lib,'SearchDepth',8,'Type','block');
hit=r(contains(r,mustContain) & contains(r,mustContain2));
assert(~isempty(hit),'block with [%s]&[%s] not found in %s',mustContain,mustContain2,lib);
s=hit{1};
end

function cleanupArtifacts(base,genlib,name,outdir)
for b={genlib,name}; if bdIsLoaded(b{1}); close_system(b{1},0); end; end
files={[base '.mat'],[genlib '.slx'],[name '.slx'],fullfile(outdir,[name '.slx'])};
for f=files; if exist(f{1},'file'); try, delete(f{1}); catch, end; end; end
dirs={['+' base],['+' base 'LumpingAdapters'],['+' genlib 'LumpingAdapters']};
for d=dirs; if exist(d{1},'dir'); try, rmdir(d{1},'s'); catch, end; end; end
end
