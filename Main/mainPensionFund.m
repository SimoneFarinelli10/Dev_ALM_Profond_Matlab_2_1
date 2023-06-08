%This is the main program generating a pension fund ALM analysis
clear;
close all

display('Starting ALM analysis...');

root = '/Users/simone/Dropbox/BackUp_Profond/ALM_Tools/Dev_ALM_Profond_Matlab_2_1/';

excelInputs = {...
                {'', 'AssetScenario.xlsx', 'LiabilityScenario.xlsx', 'CorrelationScenario.xlsx', 'ControlScenario.xlsx'},...
%                {'Sensitivity_Profond_ALM_2019_TZ_275_200_AGHZ_100\', 'AssetScenario_Profond.xlsx', 'LiabilityScenario_Profond.xlsx', 'CorrelationScenario_Profond.xlsx', 'ControlScenario_Profond.xlsx'},...
              
             };



for j=1:length(excelInputs)
    tic;
    inputFolder = [root, excelInputs{j}{1}, 'Input/'];
    outputFolder = [root, excelInputs{j}{1}, 'Output/'];

    disp(sprintf('inputFolder is %s', inputFolder));
    disp(sprintf('outputFolder is %s', outputFolder));

    assetScenario = [inputFolder, excelInputs{j}{2}]; 
    liabilityScenario = [inputFolder, excelInputs{j}{3}];
    correlationScenario = [inputFolder, excelInputs{j}{4}];
    controlScenario = [inputFolder, excelInputs{j}{5}];
    
    display('Input files...');
    disp(assetScenario);
    disp(liabilityScenario);
    disp(correlationScenario);
    disp(controlScenario);

    display('Preparing names of output files...');
    varInputsStruct0 = prepareOptionsCalcPensionFund(outputFolder);

    display('Starting numerical computation...')
    pfp=calcPensionFund(assetScenario,...
                        liabilityScenario,...
                        correlationScenario,...
                        controlScenario,...
                        varInputsStruct0);
    toc;

    %cd([root, excelInputs{1}{1}]) 
    %system(['pdflatex ', [root, excelInputs{1}{1}, 'Profond_ALM_Report.tex']])

end