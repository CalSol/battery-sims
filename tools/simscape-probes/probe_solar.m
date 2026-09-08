load_system('ee_lib');
function ok = trycase(tag, plusPort, minusPort, irrPort)
  ok=false; nm=['sc_' tag];
  if bdIsLoaded(nm); close_system(nm,0); end
  new_system(nm); load_system(nm);
  add_block('ee_lib/Sources/Solar Cell',[nm '/PV'],'Position',[200 200 240 260]);
  add_block('fl_lib/Electrical/Electrical Elements/Resistor',[nm '/R'],'Position',[320 210 360 240],'R','2');
  add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference',[nm '/G'],'Position',[200 320 230 350]);
  add_block('nesl_utility/Solver Configuration',[nm '/S'],'Position',[120 320 160 350]);
  ph=get_param([nm '/PV'],'PortHandles');
  allc=[ph.LConn ph.RConn];   % order: LConn1,LConn2,RConn1 -> but grouped
  % Build an index list matching Type order LConn1,LConn2,RConn1
  P=struct('LConn1',ph.LConn(1),'LConn2',ph.LConn(2),'RConn1',ph.RConn(1));
  rr=get_param([nm '/R'],'PortHandles'); gg=get_param([nm '/G'],'PortHandles');
  sc=get_param([nm '/S'],'PortHandles');
  gp=gg.LConn; if isempty(gp); gp=gg.RConn; end
  sp=sc.RConn; if isempty(sp); sp=sc.LConn; end
  try
    add_line(nm, P.(plusPort), rr.LConn(1),'autorouting','on');
    add_line(nm, P.(minusPort), rr.RConn(1),'autorouting','on');
    add_line(nm, P.(minusPort), gp(1),'autorouting','on');
    add_line(nm, P.(minusPort), sp(1),'autorouting','on');
    if ~isempty(irrPort)
      add_block('simulink/Sources/Constant',[nm '/Ir'],'Position',[80 190 110 210],'Value','1000');
      add_block('nesl_utility/Simulink-PS Converter',[nm '/C'],'Position',[130 190 155 210]);
      add_line(nm,'Ir/1','C/1');
      cc=get_param([nm '/C'],'PortHandles');
      add_line(nm, cc.RConn(1), P.(irrPort),'autorouting','on');
    end
    set_param(nm,'StopTime','0.01','Solver','ode23t');
    sim(nm); fprintf('[%s] SIM OK (+=%s -=%s irr=%s)\n', tag, plusPort, minusPort, irrPort);
    ok=true;
  catch ME
    fprintf('[%s] FAIL: %s\n', tag, ME.message);
  end
  close_system(nm,0);
end
trycase('A','LConn1','RConn1','LConn2');   % +=TL, -=Bot, irr=TR
trycase('B','LConn2','RConn1','LConn1');   % +=TR, -=Bot, irr=TL
trycase('C','LConn1','RConn1','');         % no irr port (fixed irradiance)
