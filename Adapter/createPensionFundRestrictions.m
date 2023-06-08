function result = createPensionFundRestrictions(assetScenarioFileName, restrType)
%createPensionFundLiabilities creates an object of class PensionFundLiabilities with inputs of the xlsx file assetScenarioFileName
%and the string restrType

num = xlsread(assetScenarioFileName, restrType, 'D3');
numNonEl = xlsread(assetScenarioFileName, restrType, 'D6');
[~, assCat0, ~] = xlsread(assetScenarioFileName, restrType, strcat('A11:A',num2str(11 + num -1)));
[matDat0, ~,  ~] = xlsread(assetScenarioFileName, restrType, strcat('B11:', char(double('D')+ numNonEl -1), num2str(11 + num)));  

result = PensionFundRestrictions(assCat0, matDat0);

end

