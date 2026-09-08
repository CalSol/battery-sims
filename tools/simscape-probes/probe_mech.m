load_system('fl_lib');
r=find_system('fl_lib','SearchDepth',6,'Type','block');
for kw={'Torque Source','Inertia','Damper','Rotational Reference','Angular Velocity'}
  h=r(contains(r,kw{1}));
  fprintf('=== %s ===\n', kw{1}); disp(h(1:min(end,4)));
end
% DC Motor mechanical port identity: check port labels
load_system('ee_lib'); new_system('dm'); load_system('dm');
add_block('ee_lib/Electromechanical/Brushed Motors/DC Motor','dm/M','Position',[100 100 160 160]);
% port line handles + names via LibraryLinkName? use get_param port labels
try
  pl=get_param('dm/M','CompiledPortDataTypes'); catch; end
% Try to read the block's port names from its mask/appearance
pc=get_param('dm/M','PortConnectivity');
for i=1:numel(pc); fprintf('DCMotor port %d Type=%s pos=[%d %d]\n', i, pc(i).Type, pc(i).Position(1), pc(i).Position(2)); end
close_system('dm',0);
