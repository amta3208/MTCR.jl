function pressure = readPostNormalShockPressure(caseFolderPath)
% Read the post normal shock pressure from normal.dat located in an
% inputted caseFolderPath

normalShockPath = makePath(caseFolderPath,'output','normal.dat');

if isfile(normalShockPath)
  fh = fopen(normalShockPath);
  while ~feof(fh)
    line = fgetl(fh);
    if contains(line,'After normal shock')
        line = fgetl(fh);
        pressure = sscanf(line,' pres  =  %e Torr');
    end
  end
else
  pressure = nan;
end
  
  
end

