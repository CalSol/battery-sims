v=ver; names={v.Name};
fprintf('Simscape Electrical present: %d\n', any(strcmp(names,'Simscape Electrical')));
fprintf('license checkout Simscape_Electrical: %d\n', license('test','Simscape_Electrical'));
% try to load the electrical library
libs={'ee_lib','elec_lib','spe_lib'};
for L=libs
  try, load_system(L{1}); fprintf('loaded %s OK\n',L{1}); catch ME; fprintf('load %s: %s\n',L{1},ME.message); end
end
% search for solar and motor blocks across loaded electrical libs
for L=libs
  if bdIsLoaded(L{1})
    r=find_system(L{1},'SearchDepth',8,'Type','block');
    sc=r(contains(lower(r),'solar')|contains(lower(r),'pv array'));
    mo=r(contains(lower(r),'dc motor')|contains(lower(r),'motor & drive')|contains(lower(r),'motor and drive'));
    fprintf('=== %s: solar blocks ===\n',L{1}); disp(sc(1:min(end,8)));
    fprintf('=== %s: motor blocks ===\n',L{1}); disp(mo(1:min(end,8)));
  end
end
