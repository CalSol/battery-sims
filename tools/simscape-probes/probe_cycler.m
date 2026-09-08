load_system('batt_lib');
% Cycler ports + dialog params
new_system('pc'); load_system('pc');
add_block('batt_lib/Cyclers/Cycler','pc/Cy','Position',[100 100 160 160]);
pc=get_param('pc/Cy','PortConnectivity');
fprintf('--- Cycler ports ---\n');
for i=1:numel(pc); fprintf(' port %d Type=%s\n', i, pc(i).Type); end
try
  dp=get_param('pc/Cy','DialogParameters'); fn=fieldnames(dp);
  fprintf('--- Cycler dialog params (%d) ---\n', numel(fn));
  disp(fn(1:min(end,30)));
catch ME; fprintf('dp err %s\n',ME.message); end
% mask names
try; mn=get_param('pc/Cy','MaskNames'); fprintf('MaskNames:\n'); disp(mn(1:min(end,30))); catch; end
close_system('pc',0);
% search cooling plate
r=find_system('batt_lib','SearchDepth',6,'Type','block');
cp=r(contains(lower(r),'cooling')|contains(lower(r),'plate'));
fprintf('--- cooling blocks in batt_lib ---\n'); disp(cp);
