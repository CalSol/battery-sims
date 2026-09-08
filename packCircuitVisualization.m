%% ========================================================================
%  packCircuitVisualization.m
%  Simulink / Simscape equivalent-electrical models that VISUALIZE the
%  36s9p pack circuit structure:
%
%    (1) ECM_cell_circuit - 2nd-order Thevenin equivalent circuit of one cell
%                           (OCV source, R0, R1||C1, R2||C2) as a real
%                           Simscape electrical schematic.
%    (2) pack_topology    - series/parallel architecture: a 2s9p module (18
%                           cells); 18 modules chained as two blocks of 9
%                           (with an air gap) => 36s9p, 129.6 V.
%    (3) pack_layout      - physical 18x18 cell arrangement: two 9x18 blocks
%                           with the central air gap, cells grouped into
%                           9x2 (2s9p) modules.
%
%  Diagrams exported as PNGs to ./sim_outputs.
%  Requires: Simulink, Simscape.
%  ========================================================================

function packCircuitVisualization()
close all; clc;
outdir = fullfile(pwd,'sim_outputs');
if ~exist(outdir,'dir'); mkdir(outdir); end
fprintf('=== Building 36s9p equivalent-circuit visualizations ===\n');
try, buildECMcell(outdir);      catch ME; fprintf(2,'ECM cell: %s\n',ME.message); end
try, buildPackTopology(outdir); catch ME; fprintf(2,'Topology: %s\n',ME.message); end
try, buildPackLayout(outdir);   catch ME; fprintf(2,'Layout: %s\n',ME.message); end
fprintf('=== Done. Visualizations written to %s ===\n', outdir);
end

%% ------------------------------------------------------------------------
function buildECMcell(outdir)
% Real Simscape schematic of the 2nd-order Thevenin cell model.
%   Controlled Current Source ports (verified): LConn1=+, RConn1=PS command,
%   RConn2=- ;  Simulink-PS Converter: Inport1=Simulink, RConn1=PS out.
name='ECM_cell_circuit';
if bdIsLoaded(name); close_system(name,0); end
new_system(name); load_system(name);

add_block('fl_lib/Electrical/Electrical Sources/DC Voltage Source',[name '/OCV'], ...
          'Position',[40 180 80 240],'Orientation','up','v0','3.6');
add_block('fl_lib/Electrical/Electrical Elements/Resistor',[name '/R0'], ...
          'Position',[140 120 200 150],'R','0.012');
add_block('fl_lib/Electrical/Electrical Elements/Resistor',[name '/R1'], ...
          'Position',[260 90 320 120],'R','0.008');
add_block('fl_lib/Electrical/Electrical Elements/Capacitor',[name '/C1'], ...
          'Position',[260 170 320 200],'c','1500');
add_block('fl_lib/Electrical/Electrical Elements/Resistor',[name '/R2'], ...
          'Position',[380 90 440 120],'R','0.006');
add_block('fl_lib/Electrical/Electrical Elements/Capacitor',[name '/C2'], ...
          'Position',[380 170 440 200],'c','20000');
add_block('fl_lib/Electrical/Electrical Sources/Controlled Current Source',[name '/I_load'], ...
          'Position',[520 150 560 210],'Orientation','down');
add_block('simulink/Sources/From Workspace',[name '/Iload'],'Position',[610 160 680 200]);
add_block('nesl_utility/Simulink-PS Converter',[name '/S2PS'],'Position',[580 165 605 195]);
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference',[name '/Gnd'], ...
          'Position',[60 300 90 330]);
add_block('nesl_utility/Solver Configuration',[name '/Solver'],'Position',[150 300 200 330]);

h=@(b) get_param([name '/' b],'PortHandles');
oc=h('OCV'); r0=h('R0'); r1=h('R1'); c1=h('C1'); r2=h('R2'); c2=h('C2');
il=h('I_load'); s2=h('S2PS'); gd=h('Gnd'); so=h('Solver');
gp=gd.LConn; if isempty(gp); gp=gd.RConn; end
sp=so.RConn; if isempty(sp); sp=so.LConn; end
ocp=oc.LConn; ocn=oc.RConn;                 % OCV up: LConn=+ (top), RConn=- (bottom)

add_line(name, ocp(1), r0.LConn(1),'autorouting','on');   % OCV+ -> R0
add_line(name, r0.RConn(1), r1.LConn(1),'autorouting','on'); % node A
add_line(name, r0.RConn(1), c1.LConn(1),'autorouting','on');
add_line(name, r1.RConn(1), r2.LConn(1),'autorouting','on'); % node B
add_line(name, r1.RConn(1), c1.RConn(1),'autorouting','on');
add_line(name, r1.RConn(1), c2.LConn(1),'autorouting','on');
add_line(name, r2.RConn(1), c2.RConn(1),'autorouting','on'); % node C (+term)
add_line(name, r2.RConn(1), il.LConn(1),'autorouting','on'); % -> load +
add_line(name, il.RConn(2), gp(1),'autorouting','on');       % load - -> ground
add_line(name, ocn(1), gp(1),'autorouting','on');            % OCV - -> ground
add_line(name, ocn(1), sp(1),'autorouting','on');            % solver on ref node
add_line(name,'Iload/1','S2PS/1');
add_line(name, s2.RConn(1), il.RConn(1),'autorouting','on'); % PS command -> RConn1

anno(name,'2nd-order Thevenin cell ECM  (OCV + R0 + R1||C1 + R2||C2)',[120 30],13,'b');
anno(name,'OCV(SOC)',[25 250],9,'k');
anno(name,'R0 ohmic',[145 100],9,'k');
anno(name,'R1||C1 charge-transfer',[250 210],9,'k');
anno(name,'R2||C2 diffusion',[380 210],9,'k');
anno(name,'Cell + terminal',[470 120],9,'r');
anno(name,'motor current draw',[600 205],9,'gray');

f=fullfile(outdir,'fig6_ECM_cell_circuit.png');
print(['-s' name],'-dpng','-r150',f);
save_system(name, fullfile(outdir,[name '.slx']));
close_system(name,0);
fprintf('  ECM cell schematic -> fig6_ECM_cell_circuit.png\n');
end

%% ------------------------------------------------------------------------
function buildPackTopology(outdir)
% Series/parallel architecture: each source = one 2s9p module (18 cells,
% 7.2 V). 18 modules = two blocks of 9 (air gap between) chained in series.
name='pack_topology';
if bdIsLoaded(name); close_system(name,0); end
new_system(name); load_system(name);

nPerRow=9; x0=70; dx=95; yTop=110; yBot=290;   % two rows (blocks) + gap
modH=55; modW=44;
mods=cell(1,18);
for m=1:18
    if m<=nPerRow; row=1; k=m;    y=yTop;
    else;          row=2; k=m-nPerRow; y=yBot; end
    x = x0 + (k-1)*dx;
    bn=sprintf('%s/M%d',name,m);
    add_block('fl_lib/Electrical/Electrical Sources/DC Voltage Source', bn, ...
              'Position',[x y x+modW y+modH],'Orientation','up','v0','7.2');
    mods{m}=get_param(bn,'PortHandles');   % up: LConn=+ top, RConn=- bottom
end
% series chain: block1 left->right (M1..M9), then down to block2, right->left
% (M10..M18) so the electrical path snakes; connect - of one to + of next.
order=[1:9, 10:18];
for i=1:numel(order)-1
    a=mods{order(i)}; b=mods{order(i+1)};
    add_line(name, a.RConn(1), b.LConn(1),'autorouting','on');  % (-) -> (+)
end

add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference',[name '/Gnd'], ...
          'Position',[x0-45 yBot+90 x0-15 yBot+120]);
gd=get_param([name '/Gnd'],'PortHandles'); gp=gd.LConn; if isempty(gp); gp=gd.RConn; end
add_line(name, mods{18}.RConn(1), gp(1),'autorouting','on');    % pack - to ref

anno(name,'36s9p pack architecture',[x0 30],14,'b');
anno(name,'each source = one 2s9p module (18 cells, 7.2 V)',[x0 70],10,'dg');
anno(name,'BLOCK 1: 9 modules  (rows going back)',[x0 yTop-30],9,'gray');
anno(name,'---- air gap between blocks ----',[x0+3*dx yBot-40],10,'b');
anno(name,'BLOCK 2: 9 modules',[x0 yBot-18],9,'gray');
anno(name,'18 modules x 2s = 36 series ; 9p  ->  129.6 V, 40.5 Ah, 324 cells', ...
     [x0 yBot+80],10,'r');
anno(name,'Pack +',[x0-30 yTop-5],10,'r');

f=fullfile(outdir,'fig7_pack_topology.png');
print(['-s' name],'-dpng','-r150',f);
save_system(name, fullfile(outdir,[name '.slx']));
close_system(name,0);
fprintf('  Pack topology schematic -> fig7_pack_topology.png\n');
end

%% ------------------------------------------------------------------------
function buildPackLayout(outdir)
% Physical 18x18 cell arrangement: two 9x18 blocks + central air gap, cells
% grouped into 9-row x 2-col (2s9p) modules (colour-coded).
D=0.0210;                       % 21700 diameter (m), for aspect only
nR=18; nC=18; gapRow=9; pitch=1.0; gapSize=1.4;
f=figure('Color','w','Position',[100 100 980 900],'Visible','off');
try, theme(f,'light'); catch, end
ax=axes(f); hold(ax,'on'); axis(ax,'equal'); set(ax,'YDir','reverse');

modCols=2;                       % each module spans 2 columns (2s)
cmap=[0.80 0.88 1.0; 0.62 0.78 0.98];   % alternating module shades
for r=1:nR
    % vertical position with a gap after gapRow
    yr = r; if r>gapRow, yr = r + gapSize; end
    blockRow = (r<=gapRow);      %#ok<NASGU>
    for cc=1:nC
        modIdx = ceil(cc/modCols);            % 1..9 module column within block
        shade  = cmap(mod(modIdx,2)+1,:);
        x=cc*pitch; y=yr*pitch;
        rectangle(ax,'Position',[x-0.45 y-0.45 0.9 0.9],'Curvature',[1 1], ...
                  'FaceColor',shade,'EdgeColor',[0.35 0.35 0.4],'LineWidth',0.4);
    end
end
% module boundary boxes (2-col wide, 9-row tall) for both blocks
for blk=0:1
    r1 = blk*gapRow + 1; r2 = r1+gapRow-1;
    y1 = r1; if r1>gapRow; y1=r1+gapSize; end
    y2 = r2; if r2>gapRow; y2=r2+gapSize; end
    for mc=1:9
        c1=(mc-1)*modCols+1;
        rectangle(ax,'Position',[c1*pitch-0.5 y1-0.5 modCols*pitch (y2-y1)+1], ...
                  'EdgeColor',[0.15 0.15 0.5],'LineWidth',1.4);
    end
end
% gap annotation
ytext=(gapRow+gapSize/2+0.5);
text(ax,nC/2+0.5, ytext, 'central air gap  (induced-draft exhaust  \leftarrow air drawn inward)', ...
     'HorizontalAlignment','center','Color',[0.1 0.1 0.6],'FontSize',10,'FontWeight','bold');
annulusArrow(ax,nC/2+0.5, 2, ytext-0.7);          % inward arrows block 1
annulusArrow(ax,nC/2+0.5, nR+gapSize-1, ytext+0.7);% inward arrows block 2

title(ax,{'36s9p Pack - Physical 18\times18 Cell Layout', ...
      'two 9\times18 blocks + central air gap; each 9\times2 outline = one 2s9p module'}, ...
      'FontSize',13,'FontWeight','bold');
xlabel(ax,'cell column (1-18)   \rightarrow   2 columns = 2s of a module');
ylabel(ax,'cell row (1-18, 9 per block)');
xlim(ax,[0.2 nC+0.8]); ylim(ax,[0.2 nR+gapSize+0.8]);
set(ax,'XTick',1:nC,'YTick',[]);
f2=fullfile(outdir,'fig8_pack_layout.png');
exportgraphics(f,f2,'Resolution',150);
close(f);
fprintf('  Physical 18x18 layout -> fig8_pack_layout.png\n');
end

function annulusArrow(ax,x,yTail,yHead)
% small vertical arrow indicating inward air flow
plot(ax,[x x],[yTail yHead],'-','Color',[0.1 0.45 0.85],'LineWidth',1.2);
dir=sign(yHead-yTail);
plot(ax,[x-0.25 x x+0.25],[yHead-0.3*dir yHead yHead-0.3*dir],'-', ...
     'Color',[0.1 0.45 0.85],'LineWidth',1.2);
end

%% ------------------------------------------------------------------------
function anno(model, txt, pos, fs, col)
% Simulink annotations take a color-NAME string (not an RGB triple).
a = Simulink.Annotation(model, txt);
a.Position = pos; a.FontSize = fs;
map=containers.Map({'b','k','r','g','gray','dg'}, ...
    {'blue','black','red','green','gray','darkGreen'});
c=char(col); if isKey(map,c); c=map(c); end
try, a.ForegroundColor=c; catch, a.ForegroundColor='black'; end
end
