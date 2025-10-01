function caseData = postProcessMTCR(casePath, caseSaveName, debug, overwrite)
% Postprocess an MTCR simulation. Will save the output to a .mat file and,
% if the user requests an output, the postprocessed simulation results will
% be returned as well.

  %% Set debug and overwrite variable if unset
  if (nargin < 3) 
    debug = true;
  end
  if (nargin < 4)
    overwrite = false;
  end

  %% Make paths
  inputPath = sprintf('%s/input',casePath);
  outputPath = sprintf('%s/output',casePath);
  statesPath = sprintf('%s/output/states',casePath);
  sourcePath = sprintf('%s/output/sources',casePath);
  
  %% Read setup data
  
  fprintf('Reading prob_setup.inp... \n');
  prob_setup = readProbSetupMTCR(inputPath);
  databasePath = prob_setup.DATABASE_PATH;
  
  fprintf('Reading species.dat... \n');
  species = readSpeciesMTCR(databasePath);
  
  %% Read output data
  
  % Initializations
  units = struct();
  labels = struct();
  
  fprintf('Reading excited state data...\n');
  [levels,species] = readLevelsMTCR(statesPath,species);
  
  % Read results
  fprintf('Reading result files...\n');
  [result, units, labels] = readResultMTCR(units,labels,outputPath,species,levels,prob_setup);
  result = trimStructMTCR(result,debug);    % Remove extra from the end of files to make them all the same length
  
  % Read sources
  fprintf('Reading source files...\n');
  [source, labels, units] = readSourceMTCR(units,labels,sourcePath,species,levels,prob_setup);
  source = trimStructMTCR(source,debug);    % Remove extra from the end of files to make them all the same length
  
  %% Save run file 
  % Save the output to a .mat file in the case directory
  saveOutputsMTCR(casePath, caseSaveName, result, source, units, labels, levels, species, prob_setup, overwrite);

  % Assign output if one is requested
  if (nargout == 1)
    caseData = struct('result',result,'source',source,'units',units,'labels',labels,'species',species,'prob_setup',prob_setup);
  end
end
