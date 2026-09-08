%% ========================================================================
%  thermalRunawayStudy_36s9p.m
%  Single-cell thermal-runaway FAILURE ANALYSIS for the 36s9p 21700 pack
%  (physical 18x18 grid: two 9x18 blocks separated by the central air gap).
%
%  Objective:
%    One cell in block 1 is driven into thermal runaway (TR) by an internal-
%    short fault.  Determine (a) whether TR propagates to neighbours and
%    whether the inter-block air gap acts as a firebreak, (b) the minimum
%    FORCED-AIR cooling required to arrest within-block propagation, and
%    (c) how those requirements relate to realistic BMS temperature-sensing
%    thresholds and the time available for the BMS to detect and react.
%
%  Physics:
%    * Lumped finite-volume thermal network of the 18x18 cell grid with the
%      central gap (no conduction OR radiation across it) -> firebreak.
%    * TR heat release: single-equation Arrhenius decomposition with reactant
%      depletion (Hatchard/NREL-style), ~150 C onset, ~45 kJ/cell (peak>600C).
%    * Neighbour ignition when a cell crosses the onset temperature.
%    * Cooling = negative-pressure forced air (convection coefficient swept).
%
%  Requires: MATLAB (base). Consistent with batterySimFramework_36s9p.m.
%  ========================================================================

function thermalRunawayStudy_36s9p()

close all; clc;
outdir = fullfile(pwd,'sim_outputs');
if ~exist(outdir,'dir'); mkdir(outdir); end

P = studyParams();
fprintf('==================================================================\n');
fprintf(' 36s9p Thermal-Runaway Failure Analysis (18x18, central air gap)\n');
fprintf('==================================================================\n');
fprintf(' Trigger cell: block-1 (r=%d,c=%d) | onset %.0f C | TR energy %.0f kJ/cell\n', ...
        P.trig(1),P.trig(2),P.T_onset_C,P.E_TR/1e3);

%% ---------------- PART A: baseline propagation at nominal cooling -------
fprintf('\n[A] Baseline single-cell TR at nominal air cooling (h=%.0f W/m^2K)...\n', P.h_nom);
base = runTR(P, P.h_nom);
fprintf('    Trigger peak: %.0f C at t=%.0fs\n', base.Ttrig_peak, base.t_trig_peak);
fprintf('    Hottest neighbour peak: %.0f C\n', base.Tnb_peak);
if base.propagated
    fprintf('    >> PROPAGATION: %d cells exceeded onset (within-block cascade).\n', base.n_ignited);
else
    fprintf('    >> CONTAINED: no neighbour reached onset.\n');
end
fprintf('    Block-2 (across gap) peak: %.0f C  -> gap firebreak %s\n', ...
        base.Tblk2_peak, ternary(base.Tblk2_peak<P.T_onset_C,'HOLDS','BREACHED'));

%% ---------------- PART B: cooling sweep -> minimum to arrest ------------
fprintf('\n[B] Sweeping forced-air coefficient to find arrest threshold...\n');
hsweep = P.hsweep;
Tnb = zeros(size(hsweep)); prop = false(size(hsweep)); nign = zeros(size(hsweep));
Qrej = zeros(size(hsweep));
for i=1:numel(hsweep)
    r = runTR(P, hsweep(i));
    Tnb(i)=r.Tnb_peak; prop(i)=r.propagated; nign(i)=r.n_ignited; Qrej(i)=r.Qcool_peak;
    fprintf('    h=%3.0f W/m^2K : neighbour peak %5.0f C | %s\n', hsweep(i), Tnb(i), ...
            ternary(prop(i),'PROPAGATES','contained'));
end
% minimum cooling that prevents propagation
idx_ok = find(~prop, 1, 'first');
if isempty(idx_ok)
    h_crit = NaN; fprintf('    No swept cooling level prevented propagation.\n');
else
    h_crit = hsweep(idx_ok);
    fprintf('    >> Minimum forced-air coefficient to arrest propagation: h_crit ~ %.0f W/m^2K\n', h_crit);
    fprintf('       (air must carry away %.0f W from the fault zone at that point)\n', Qrej(idx_ok));
end

%% ---------------- PART C: BMS sensing threshold analysis ----------------
fprintf('\n[C] BMS detection analysis...\n');
% (i) UNDER-COOLED case (a cooling level that still propagates): this is the
%     dangerous scenario the BMS must catch.  Use the level just below h_crit.
idx_prop = find(prop, 1, 'last');
h_uncooled = hsweep(max(idx_prop,1));
danger = runTR(P, h_uncooled);
bmsD   = bmsDetection(danger, P);
t_nb_onset = crossTime(danger.t, danger.Tnbmax, P.T_onset_C);   % propagation deadline
window = t_nb_onset - bmsD.t_first;

% (ii) DESIGN case (cooling that arrests): still detectable, safe margin.
h_design = pickDesign(hsweep, prop, P.h_nom);
det = runTR(P, h_design);
bms = bmsDetection(det, P);

fprintf('  Dangerous (under-cooled) scenario  h=%.0f W/m^2K:\n', h_uncooled);
fprintf('    BMS sensor (adjacent cell) trips: %s at t=%s s\n', bmsD.first_kind, fmt(bmsD.t_first));
fprintf('      60C warn %s s | 80C crit %s s | dT/dt %s s\n', ...
        fmt(bmsD.t_warn), fmt(bmsD.t_crit), fmt(bmsD.t_rate));
fprintf('    Neighbour would reach onset (propagation) at t=%s s\n', fmt(t_nb_onset));
fprintf('    >> BMS actionable window before cascade: %s s\n', fmt(window));
fprintf('  Design (arrested) scenario  h=%.0f W/m^2K:\n', h_design);
fprintf('    Sensor detects (%s) at t=%s s; neighbour peak stays %.0f C (<%.0f C onset)\n', ...
        bms.first_kind, fmt(bms.t_first), det.Tnb_peak, P.T_onset_C);
fprintf('    Air heat-rejection capacity required: %.0f W from the fault zone\n', det.Qcool_peak);
fprintf('  Inter-block gap: block-2 peak %.0f C -> the central air gap %s propagation.\n', ...
        det.Tblk2_peak, ternary(det.Tblk2_peak<P.T_onset_C,'blocks','does NOT block'));
fprintf('  DESIGN RULE: draw forced air giving h >= h_crit ~ %.0f W/m^2K within ~%s s of a\n', ...
        h_crit, fmt(window));
fprintf('    rate/threshold trip to prevent single-cell TR from cascading within a block.\n');

%% ---------------- Visualization ----------------------------------------
makeTRfigures(P, base, hsweep, Tnb, prop, h_crit, det, bms, danger, bmsD, ...
              t_nb_onset, h_design, h_uncooled, outdir);

save(fullfile(outdir,'thermal_runaway_results.mat'), ...
     'P','base','hsweep','Tnb','prop','h_crit','det','bms','h_design');
fprintf('\nOutputs written to: %s\n', outdir);
fprintf('==================================================================\n');
end

%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================
function P = studyParams()
% Physical 18x18 grid (two 9x18 blocks + central air gap), 21700 cells
P.nR=18; P.nC=18; P.gapRow=9;   % gap between rows 9 and 10 (firebreak)
P.trig=[5 9];                   % trigger cell (interior of block 1)
% cell thermal
P.mass=0.068; P.cp=1000; P.Cth=P.mass*P.cp;         % J/K
P.D=0.021; P.H=0.070; P.r=P.D/2;
P.A_side=2*pi*P.r*P.H;
% cell-cell conduction (line contact + air), see framework audit note
gap=0.002; kc=2.0; ka=0.026; At=0.015*P.H*P.D; Aa=P.H*P.D-At;
P.G_cc=(kc*At+ka*Aa)/gap;       % W/K
% cooling = negative-pressure forced air (convection coefficient swept)
P.Tcool_C=25; P.T_amb_C=25;
P.h_nom=60;                     % nominal forced-air coefficient (W/m^2K)
% ---- Thermal-runaway kinetics (Arrhenius, 1-equation, reactant depletion)
P.A_arr=1e13;                   % 1/s pre-exponential
P.Ea=1.30e5;                    % J/mol activation energy
P.Rgas=8.314;
P.E_TR=45e3;                    % J total exothermic energy per cell
P.T_onset_C=150;                % neighbour ignition threshold
% ---- Cell-to-cell radiative + vent-gas coupling (drives TR propagation) ---
%   During runaway the failed cell vents hot gas/ejecta onto neighbours; the
%   effective can-to-can coupling is far larger than nominal conduction.
P.eps=0.85; P.sigma=5.67e-8;    % emissivity, Stefan-Boltzmann
P.A_rad=0.5*P.A_side;           % m^2 effective facing area (dense module)
% ---- Internal-short trigger (fault that initiates TR) ----
P.P_short=250;                  % W localized short power (hard internal short)
P.t_short=150;                  % s duration of the fault
% ---- Time integration ----
P.dt=0.02; P.t_end=260;
% ---- BMS thresholds ----
P.T_warn_C=60; P.T_crit_C=80; P.dTdt_trip=2.0;   % C, C, C/s
% ---- Cooling sweep (refined near the arrest threshold) ----
P.hsweep=[20 30 40 45 50 55 60 70 85 100 150 200 300 450 600 900];
end

function S = runTR(P, h_cool)
% Simulate one-cell TR with given forced-air coefficient h_cool (W/m^2K).
nR=P.nR; nC=P.nC; N=nR*nC; dt=P.dt; gr=P.gapRow;
Cth=P.Cth;
% thermal operators: the central gap breaks BOTH conduction and radiation
% between rows gr and gr+1 (firebreak).
Lg = buildLap(nR,nC,P.G_cc,gr);
Lu = buildLap(nR,nC,1,gr);                       % unit Laplacian for radiation
kRad = P.eps*P.sigma*P.A_rad;                    % radiative coupling coeff
Gcool = h_cool*P.A_side*ones(N,1);              % W/K per cell
A = (Cth/dt)*speye(N) + Lg + spdiags(Gcool,0,N,N);
dA = decomposition(A,'lu');
Tcool = P.Tcool_C+273.15;

% indices
idx=@(r,c)(c-1)*nR+r;
it = idx(P.trig(1),P.trig(2));                  % trigger linear index
% neighbours of trigger (4-connected)
nb=[];
if P.trig(1)>1, nb(end+1)=idx(P.trig(1)-1,P.trig(2)); end
if P.trig(1)<nR,nb(end+1)=idx(P.trig(1)+1,P.trig(2)); end
if P.trig(2)>1, nb(end+1)=idx(P.trig(1),P.trig(2)-1); end
if P.trig(2)<nC,nb(end+1)=idx(P.trig(1),P.trig(2)+1); end

% state
T=(P.T_amb_C+273.15)*ones(N,1);
c=ones(N,1);                                    % reactant fraction
nt=round(P.t_end/dt)+1; t=(0:nt-1)'*dt;
Ttrig=zeros(nt,1); Tnbmax=zeros(nt,1); Tsens=zeros(nt,1);
Qcool_hist=zeros(nt,1);
Ttrig(1)=P.T_amb_C; Tnbmax(1)=P.T_amb_C; Tsens(1)=P.T_amb_C;
sens = nb(1);                                    % primary BMS sensor = adjacent cell
Tfield_peak = T;
maxT_all = (P.T_amb_C+273.15)*ones(N,1);

for k=2:nt
    tk=t(k);
    % Arrhenius reaction heat (explicit, energy-conserving reactant burn)
    kA = P.A_arr*exp(-P.Ea./(P.Rgas*T));
    dc = min(kA.*c*dt, c);                        % reactant consumed (bounded)
    c  = c - dc;
    Q_ra = P.E_TR*dc/dt;                          % W per cell (energy conserving)
    % cell-to-cell radiation (dominant propagation path at TR temperatures)
    Q_rad = -kRad*(Lu*(T.^4));                    % net radiative gain per cell
    % internal-short trigger
    Q = Q_ra + Q_rad;
    if tk<=P.t_short, Q(it)=Q(it)+P.P_short; end
    % implicit conduction+cooling step
    b = (Cth/dt)*T + Q + Gcool*Tcool;
    T = dA\b;
    maxT_all = max(maxT_all,T);
    % logs
    Ttrig(k)=T(it)-273.15;
    Tnbmax(k)=max(T(nb))-273.15;
    Tsens(k)=T(sens)-273.15;
    Qcool_hist(k)=sum(Gcool.*(T-Tcool));         % W rejected to coolant
end

S.t=t; S.Ttrig=Ttrig; S.Tnbmax=Tnbmax; S.Tsens=Tsens; S.Qcool_hist=Qcool_hist;
S.h_cool=h_cool; S.nb=nb; S.sens=sens; S.it=it;
[S.Ttrig_peak,ip]=max(Ttrig); S.t_trig_peak=t(ip);
S.Tnb_peak=max(Tnbmax);
S.n_ignited = nnz(maxT_all-273.15 >= P.T_onset_C) - 1;   % excl. trigger
S.propagated = S.Tnb_peak >= P.T_onset_C;
S.Qcool_peak = max(Qcool_hist);
S.maxT_field = reshape(maxT_all-273.15,nR,nC);
% block-2 peak (across the gap) to test the firebreak
blk2 = S.maxT_field(gr+1:end,:);
S.Tblk2_peak = max(blk2(:));
end

function b = bmsDetection(S,P)
% Detection times for warning/critical/rate thresholds on the adjacent sensor
t=S.t; Ts=S.Tsens;
b.t_warn = crossTime(t,Ts,P.T_warn_C);
b.t_crit = crossTime(t,Ts,P.T_crit_C);
dTdt=[0;diff(Ts)./diff(t)];
b.t_rate = crossTime(t,dTdt,P.dTdt_trip);
% trigger cell reaches onset (runaway becomes unstoppable ~ onset)
b.t_trig_onset = crossTime(t,S.Ttrig,P.T_onset_C);
cand=[b.t_warn b.t_crit b.t_rate];
kinds={'60C warning','80C critical','dT/dt rate'};
[b.t_first,ci]=min(cand);
if isempty(ci)||isnan(b.t_first); b.first_kind='(none)'; else; b.first_kind=kinds{ci}; end
b.lead_time = b.t_trig_onset - b.t_first;
b.dTdt=dTdt;
end

function h=pickDesign(hsweep,prop,hnom)
ok=find(~prop,1,'first');
if isempty(ok); h=hsweep(end); else; h=hsweep(ok); end
h=max(h,hnom);
end

function L=buildLap(nR,nC,g,gapRow)
% 4-neighbour grid Laplacian; if gapRow>0, no coupling between rows gapRow
% and gapRow+1 (central air gap = firebreak).
if nargin<4, gapRow=0; end
N=nR*nC; idx=@(r,c)(c-1)*nR+r; I=[];J=[];S=[];
for c=1:nC; for r=1:nR
    i=idx(r,c); nb=[];
    if r>1  && (r-1)~=gapRow, nb(end+1)=idx(r-1,c); end
    if r<nR && r~=gapRow,     nb(end+1)=idx(r+1,c); end
    if c>1,nb(end+1)=idx(r,c-1);end; if c<nC,nb(end+1)=idx(r,c+1);end
    for j=nb; I(end+1)=i;J(end+1)=j;S(end+1)=-g; end %#ok<AGROW>
    I(end+1)=i;J(end+1)=i;S(end+1)=g*numel(nb); %#ok<AGROW>
end; end
L=sparse(I,J,S,N,N);
end

function tc=crossTime(t,y,thr)
k=find(y>=thr,1,'first');
if isempty(k); tc=NaN; else; tc=t(k); end
end

function s=fmt(x); if isnan(x); s='n/a'; else; s=sprintf('%.1f',x); end; end
function y=ternary(c,a,b); if c; y=a; else; y=b; end; end

function makeTRfigures(P,base,hsweep,Tnb,prop,h_crit,det,bms,danger,bmsD, ...
                       t_nb_onset,h_design,h_uncooled,outdir)
LT=@(f) applyTheme(f);

% ---- Figure 4: TR trajectories, BMS window, cooling sweep, cascade map ----
f=figure('Color','w','Position',[60 60 1500 820],'Name','Thermal Runaway'); LT(f);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
title(tl,'Single-Cell Thermal Runaway - Failure & Air-Cooling Analysis (36s9p, 18x18)', ...
      'FontSize',14,'FontWeight','bold');

% (a) design (arrested) case: trigger runs away but neighbour is held down
nexttile; hold on; grid on;
plot(det.t,det.Ttrig,'r','LineWidth',2.2);
plot(det.t,det.Tnbmax,'Color',[0.1 0.55 0.2],'LineWidth',1.8);
plot(danger.t,danger.Tnbmax,'--','Color',[0.85 0.35 0.1],'LineWidth',1.6);
yline(P.T_onset_C,'--k','Onset 150\circC');
xlabel('Time (s)'); ylabel('Temperature (\circC)');
title(sprintf('(a) Trigger vs neighbour: adequate vs under cooling'));
legend('trigger cell (runaway)', ...
       sprintf('neighbour @ h=%.0f (arrested)',h_design), ...
       sprintf('neighbour @ h=%.0f (cascades)',h_uncooled),'Location','east');
xlim([0 P.t_end]);

% (b) BMS detection timeline in the dangerous (under-cooled) scenario
nexttile; hold on; grid on;
plot(danger.t,danger.Tsens,'Color',[0.1 0.45 0.85],'LineWidth',2);
yline(P.T_warn_C,'--','Color',[0 0.5 0],'Label','60\circC warn');
yline(P.T_crit_C,'--','Color',[0.6 0 0],'Label','80\circC crit');
mark(bmsD.t_warn,P.T_warn_C,[0 0.5 0]);
mark(bmsD.t_crit,P.T_crit_C,[0.6 0 0]);
if ~isnan(bmsD.t_rate); xline(bmsD.t_rate,':','Color',[0.2 0.2 0.8],'Label','dT/dt trip'); end
if ~isnan(t_nb_onset);  xline(t_nb_onset,'-r','Label','neighbour ignites'); end
if ~isnan(bmsD.t_first) && ~isnan(t_nb_onset)
    yl=ylim; xa=[bmsD.t_first t_nb_onset];
    patch([xa fliplr(xa)],[yl(1) yl(1) yl(2) yl(2)],[0.2 0.4 0.9], ...
          'FaceAlpha',0.08,'EdgeColor','none');
    text(mean(xa),yl(1)+0.85*diff(yl),sprintf('BMS window\n%.0f s',t_nb_onset-bmsD.t_first), ...
         'HorizontalAlignment','center','FontSize',9,'Color',[0 0 0.6]);
end
xlabel('Time (s)'); ylabel('Adjacent-sensor temp (\circC)');
title('(b) BMS detection vs propagation deadline');
xlim([0 min(t_nb_onset*1.4,P.t_end)]);

% (c) cooling sweep: neighbour peak temp vs cooling
nexttile; hold on; grid on;
okp=~prop;
plot(hsweep(okp),Tnb(okp),'o','MarkerFaceColor',[0.2 0.6 0.2],'Color',[0.2 0.6 0.2],'MarkerSize',7);
plot(hsweep(~okp),Tnb(~okp),'s','MarkerFaceColor',[0.8 0.2 0.2],'Color',[0.8 0.2 0.2],'MarkerSize',7);
plot(hsweep,Tnb,'k-','LineWidth',1.0);
yline(P.T_onset_C,'--k','Onset (propagation)');
if ~isnan(h_crit); xline(h_crit,'-.','Color',[0 0 0.7],'Label','h_{crit}'); end
set(gca,'XScale','log');
xlabel('Forced-air coefficient h (W/m^2K, log)'); ylabel('Peak neighbour temp (\circC)');
title('(c) Propagation threshold vs air cooling');
legend('contained','propagates','Location','northeast');

% (d) peak-temperature field of the UNDER-COOLED cascade, with the gap shown
gr=P.gapRow;
Fdisp=[danger.maxT_field(1:gr,:); nan(1,P.nC); danger.maxT_field(gr+1:end,:)];
nexttile; hI=imagesc(Fdisp); set(hI,'AlphaData',~isnan(Fdisp));
axis image; colormap(gca,hot); set(gca,'Color',[0.9 0.9 0.92]); cb=colorbar;
ylabel(cb,'Peak temp (\circC)');
xlabel('Cell column (1-18)'); ylabel('Cell row (1-18, gap between blocks)');
title(sprintf('(d) Cascade @ h=%.0f W/m^2K: gap firebreak (blk2 peak %.0f\\circC)', ...
      h_uncooled, danger.Tblk2_peak));
hold on; plot(P.trig(2),P.trig(1),'cp','MarkerSize',13,'MarkerFaceColor','c');
text(P.nC/2, gr+1, 'air gap (firebreak)','HorizontalAlignment','center', ...
     'FontSize',8,'Color',[0.15 0.15 0.5],'FontAngle','italic');

exportgraphics(f,fullfile(outdir,'fig4_thermal_runaway.png'),'Resolution',150);
fprintf('    Figure saved: fig4_thermal_runaway.png\n');
end

function mark(tx,ty,col)
if ~isnan(tx); plot(tx,ty,'v','MarkerSize',9,'MarkerFaceColor',col,'MarkerEdgeColor','k'); end
end
function applyTheme(f); try, theme(f,'light'); catch, end; set(f,'Color','w'); end
