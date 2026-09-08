load_system('ee_lib');
new_system('pp'); load_system('pp');
add_block('ee_lib/Sources/Solar Cell','pp/PV','Position',[200 200 240 260]);
add_block('ee_lib/Electromechanical/Brushed Motors/DC Motor','pp/M','Position',[400 200 460 260]);
for b={'PV','M'}
  pc=get_param(['pp/' b{1}],'PortConnectivity');
  ps=get_param(['pp/' b{1}],'Position');
  fprintf('=== %s  block pos [%d %d %d %d] ===\n', b{1}, ps(1),ps(2),ps(3),ps(4));
  for i=1:numel(pc)
    fprintf('  port %d Type=%-7s pos=[%d %d]\n', i, pc(i).Type, pc(i).Position(1), pc(i).Position(2));
  end
end
% Solar cell: is irradiance an input port? check block variants / ports labels
try
  ph=get_param('pp/PV','PortHandles');
  fprintf('PV: LConn=%d RConn=%d Inport=%d\n', numel(ph.LConn),numel(ph.RConn),numel(ph.Inport));
catch; end
% check Solar Cell has an irradiance port option
try
  v=get_param('pp/PV','ee_solar_cell_ipts'); fprintf('PV iports param: %s\n', v); catch; end
close_system('pp',0);
