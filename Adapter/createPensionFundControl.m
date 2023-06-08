function [projectionFlag, numSim, horProj, fLevelsShortfall, fLevelRes, alphaRes, optimizationFlag, horOpt, numOpt, retRangeOpt, fLevelOpt] = createPensionFundControl(controlScenarioFileName)
  
  projectionFlag  = xlsread(controlScenarioFileName, 'Control', 'B3');
  numSim = xlsread(controlScenarioFileName, 'Control', 'B4'); 
  horProj = xlsread(controlScenarioFileName, 'Control', 'B5');
  fLevelsShortfall = xlsread(controlScenarioFileName, 'Control', 'B6:J6');
  fLevelRes = xlsread(controlScenarioFileName, 'Control', 'B7');
  alphaRes = xlsread(controlScenarioFileName, 'Control', 'B8');

  optimizationFlag = xlsread(controlScenarioFileName, 'Control', 'B10');
  horOpt = xlsread(controlScenarioFileName, 'Control', 'B11');
  numOpt = xlsread(controlScenarioFileName, 'Control', 'B12');
  retRangeOpt = xlsread(controlScenarioFileName, 'Control', 'B13:C13');
  fLevelOpt = xlsread(controlScenarioFileName, 'Control', 'B14');
   
end

