try
    packElectricalModel_36s9p();
catch ME
    fprintf(2,'ERR: %s\n', ME.message);
    for i=1:min(numel(ME.stack),6); fprintf(2,'  at %s line %d\n', ME.stack(i).name, ME.stack(i).line); end
end