function result = createPensionFundProjection(assetScenarioFileName, liabilityScenarioFileName, correlationFileName)
%createPensionFundProjection creates an object of class PensionFundProjection with inputs of the xlsx files assetScenarioFile, 
%liabilityScenarioFile and correlationFile."

 [date, ~, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', 'B1:D1');
 num = xlsread(correlationFileName, 'Korrelationen', 'D3');
 hor = xlsread(correlationFileName, 'Korrelationen', 'D6');
[corrDat, corrText,  ~] = xlsread(correlationFileName, 'Korrelationen', GetExcelRange(10, 10 + num, 1, hor * (num + 1)));

 runningCorr = corrDat;
 corr.('Cat') = corrText(2:end,1);
 for i=1:hor
   corr.(strcat('Period_',num2str(i))) = runningCorr(1:end,1:num);
   runningCorr = runningCorr(:,num+2:end);
 end
 pfa = createPensionFundAssets(assetScenarioFileName);
 pfl = createPensionFundLiabilities(liabilityScenarioFileName);
 pfr = createPensionFundRestrictions(assetScenarioFileName, 'Restriktionen');
 pfrOpt = createPensionFundRestrictions(assetScenarioFileName, 'Optimierungsrestriktionen');
 prDate = datenum(date);
 result = PensionFundProjection(pfa, pfl, corr, pfr, pfrOpt, prDate);

end
