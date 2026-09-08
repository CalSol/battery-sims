load_system('ee_lib');
new_system('ps'); load_system('ps');
add_block('ee_lib/Sources/Solar Cell','ps/PV','Position',[100 100 140 140]);
add_block('ee_lib/Electromechanical/Brushed Motors/DC Motor','ps/M','Position',[300 100 360 150]);
for b={'PV','M'}
  pc=get_param(['ps/' b{1}],'PortConnectivity');
  fprintf('=== %s ports ===\n', b{1});
  for i=1:numel(pc); fprintf('  %d Type=%s\n', i, pc(i).Type); end
  mn=get_param(['ps/' b{1}],'MaskNames');
  fprintf('  mask(%d): %s\n', numel(mn), strjoin(mn(1:min(end,14))',', '));
end
% tiny sim: solar cell -> resistor, to confirm license checks out
add_block('fl_lib/Electrical/Electrical Elements/Resistor','ps/R','Position',[200 220 240 250],'R','5');
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference','ps/G','Position',[120 300 150 330]);
add_block('nesl_utility/Solver Configuration','ps/S','Position',[60 300 100 330]);
pv=get_param('ps/PV','PortHandles'); rr=get_param('ps/R','PortHandles');
gg=get_param('ps/G','PortHandles'); ss=get_param('ps/S','PortHandles');
gp=gg.LConn; if isempty(gp); gp=gg.RConn; end
sp=ss.RConn; if isempty(sp); sp=ss.LConn; end
% Solar Cell ports: find +/- conserving among LConn/RConn
fprintf('PV LConn=%d RConn=%d\n', numel(pv.LConn), numel(pv.RConn));
try
  add_line('ps', pv.RConn(1), rr.LConn(1),'autorouting','on');
  add_line('ps', pv.LConn(1), rr.RConn(1),'autorouting','on');
  add_line('ps', pv.LConn(1), gp(1),'autorouting','on');
  add_line('ps', pv.LConn(1), sp(1),'autorouting','on');
  set_param('ps','StopTime','0.01','Solver','ode23t');
  out=sim('ps');
  fprintf('SIM OK (Electrical license works at sim time)\n');
catch ME
  fprintf(2,'SIM/wire err: %s\n', ME.message);
end
close_system('ps',0);
