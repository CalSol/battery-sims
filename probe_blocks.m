v=ver; names={v.Name};
for w={'Simscape','Simscape Battery','Simscape Electrical','Simscape Fluids','Simscape Driveline'}
    fprintf('%-22s : %d\n', w{1}, any(strcmp(names,w{1})));
end
fprintf('--- searching for blocks ---\n');
% load battery library
try, load_system('batterylib'); catch, end
try, load_system('batt_lib'); catch, end
cands = {'Cycler','Cooling Plate','Battery (Table-Based)','Battery Module', ...
         'Battery Equivalent Circuit','Parallel Assembly'};
libs = {'batt_lib','batterylib','fl_lib','simscape'};
for L=libs
  try
    r=find_system(L{1},'SearchDepth',6,'Type','block');
    hits=r(contains(lower(r),'cycler')|contains(lower(r),'cooling plate'));
    if ~isempty(hits); fprintf('In %s:\n',L{1}); disp(hits(1:min(end,12))); end
  catch ME
    fprintf('lib %s: %s\n', L{1}, ME.message);
  end
end
% thermal liquid pump / source blocks (foundation)
try
  rp=find_system('fl_lib','SearchDepth',6,'Type','block');
  pumps=rp(contains(lower(rp),'pump')|contains(lower(rp),'mass flow'));
  fprintf('--- fl_lib pump/flow blocks ---\n'); disp(pumps(1:min(end,15)));
catch ME; fprintf('fl_lib pump search: %s\n',ME.message); end
