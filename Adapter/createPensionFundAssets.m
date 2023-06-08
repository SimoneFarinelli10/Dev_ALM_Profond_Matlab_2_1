function result = createPensionFundAssets(assetScenarioFileName)
%createPensionFundAssets creates an object of class PensionFundAssets with inputs of the xlsx file assetScenarioFile.";

hor = xlsread(assetScenarioFileName, 'Returns', 'D3');
num = xlsread(assetScenarioFileName, 'Returns', 'D6');
numTest1 = xlsread(assetScenarioFileName, 'Anlageseite', 'C6');
numTest2 =  xlsread(assetScenarioFileName, 'Alternative Strategien', 'B6');

[~, assCat, ~] = xlsread(assetScenarioFileName, 'Returns', strcat('B10:B',num2str(10 + 4 * num -1)));
assCat = assCat(not(cellfun(@(a)strcmp(a,''), assCat)));

[~, testList1, ~] = xlsread(assetScenarioFileName, 'Anlageseite', strcat('B9:B',num2str(9 + numTest1 -1)));
if not(isequal(testList1, assCat))
    assCat
    testList1
    msgID = 'createPensionFundAssets:badargException';
    msg = 'Incompatibility between asset list in return and value sheets.';
    throw(MException(msgID,msg));
end

[~, testList2, ~] = xlsread(assetScenarioFileName, 'Alternative Strategien', strcat('A11:A',num2str(10 + numTest2)));
if not(isequal(testList2, assCat))
    assCat
    testList2
    msgID = 'createPensionFundAssets:badargException';
    msg = 'Incompatibility between asset list in return and alternative strategy sheets.';
    throw(MException(msgID,msg));
end

numAlt = xlsread(assetScenarioFileName, 'Alternative Strategien', 'B8');

assDat = xlsread(assetScenarioFileName, 'Returns', strcat('D10:', char(double('D')+ hor -1), num2str(10 + 4 * num -1)));
markV0 =  xlsread(assetScenarioFileName, 'Anlageseite', strcat('C9:C',  num2str(9 + num -1)));
[~, ~, markVEvolutionRaw] = xlsread(assetScenarioFileName, 'Anlageseite', strcat('D9:', char(double('D')+ hor - 1), num2str(9 + num - 1)));
posNum = cellfun(@isnumeric, markVEvolutionRaw);
markVEvolution = NaN(num,hor);
mat = cell2mat(markVEvolutionRaw(posNum));
markVEvolution(posNum) = mat;

[~, ~, altStratRaw] = xlsread(assetScenarioFileName, 'Alternative Strategien',  GetExcelRange(10, 11 + num , 1,  numAlt * (2 + hor)));
[~, ~, altStratNames] = xlsread(assetScenarioFileName, 'Alternative Strategien',  GetExcelRange(10, 10 , 1,  numAlt * (2 + hor)));
posNum = cellfun(@isnumeric, altStratRaw);
altStratMat = NaN(num + 2, numAlt * (2 + hor)+1);
mat = cell2mat(altStratRaw(posNum));
altStratMat(posNum) = mat;
runningAltStratMat = altStratMat;

altStratNames = altStratNames(mod([1:length(altStratNames)], hor+2)==2);


for i=1:numAlt
  altStrat.(altStratNames{i}) = runningAltStratMat(2:end,2:(hor+2));
  runningAltStratMat = runningAltStratMat(:,hor+3:end);
end
[~, famCat, ~] = xlsread(assetScenarioFileName, 'Anlageseite', strcat('A9:A',num2str(9 + num -1)));
posCat = find(not(cellfun(@(a)strcmp(a,''), famCat)));
%famCat = famCat(not(cellfun(@(a)strcmp(a,''), famCat)));
secCat = assCat(posCat(2):(posCat(3)-1));
costs = xlsread(assetScenarioFileName, 'Returns', strcat('D', num2str(10 + 4 * num + 2), ':', char(double('D')+ hor - 1), num2str(10 + 4 * num + 2))); 
realEstateCosts = xlsread(assetScenarioFileName, 'Anlageseite', strcat('D', num2str(9 + num),':', char(double('D')+ hor - 1), num2str(9 + num))); 

result = PensionFundAssets(assCat, assDat, markV0, markVEvolution, altStrat, secCat, costs, realEstateCosts);

end

