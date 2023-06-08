function result = createPensionFundLiabilities(liabilityScenarioFileName)
%createPensionFundLiabilities creates an object of class PensionFundLiabilities with inputs of the xlsx file liabilityScenarioFileNameName.';

   hor = xlsread(liabilityScenarioFileName, 'Liabilities', 'D3');
   runOff = xlsread(liabilityScenarioFileName, 'Run Off', 'D4');
   numILCat = xlsread(liabilityScenarioFileName, 'Liabilities', 'D6');
   numNILCat = xlsread(liabilityScenarioFileName, 'Liabilities', 'D7');
   numCIFCat = xlsread(liabilityScenarioFileName, 'Liabilities', 'D8');
   numCOFCat = xlsread(liabilityScenarioFileName, 'Liabilities', 'D9');
   numPsLiab = xlsread(liabilityScenarioFileName, 'Returns', 'D6');
   numQuant = xlsread(liabilityScenarioFileName, 'Quantiles', 'D4');
   numTest1 = xlsread(liabilityScenarioFileName, 'Run Off', 'D7');
   numTest2 = xlsread(liabilityScenarioFileName, 'Run Off', 'D8');
   numTest3 = xlsread(liabilityScenarioFileName, 'Returns', 'D6');
   
   [~, iLCat, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('A13:A',num2str(13 + numILCat -1)));
   [~, iCashIFCat, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('A', num2str(13 + numILCat + 3), ':A', num2str(13 + numILCat + 2 + numCIFCat)));
   [~, iCashOFCat, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('A', num2str(13 + numILCat + numCIFCat + 5), ':A', num2str(13 + numILCat + numCIFCat + numCOFCat + 4)));
   [~, nILCat, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('A', num2str(13 + numILCat + numCIFCat + numCOFCat + 7), ':A', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 6)));
   [ilDat, ~,  ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B13:', char(double('B')+ hor), num2str(13 + numILCat -1)));
   [iCashIFDat, ~, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B',num2str(13 + numILCat + 3), ':', char(double('B')+ hor -1), num2str(13 + numILCat + 2 + numCIFCat)));
   [iCashOFDat, ~,  ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B',num2str(13 + numILCat + numCIFCat + 5), ':', char(double('B')+ hor -1), num2str(13 + numILCat + numCIFCat + numCOFCat + 4)));
   [nilDat, ~, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B',num2str(13 + numILCat + numCIFCat + numCOFCat + 7), ':C', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 6)));
   [~, cfCheck, ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('C', num2str(13 + numILCat + numCIFCat + numCOFCat + 7), ':C', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 6)));
   pos = ismember(cfCheck, 'ja');
   nILCfRelCat = nILCat(pos);
   [~, psCat, ~]= xlsread(liabilityScenarioFileName, 'Returns', strcat('B10:B', num2str(10 +  4 * numPsLiab -1)));
   psCat = psCat(not(cellfun(@(a)strcmp(a,''), psCat)));
   [psDat, ~,  ~] = xlsread(liabilityScenarioFileName, 'Returns', strcat('D10:', char(double('D')+ hor -1), num2str(10 + 4 * numPsLiab -1)));
   [vfRes, ~,  ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 8),':B', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 8)));
   [mvDev, ~,  ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B', num2str(13 + numILCat + numCIFCat + numCOFCat + numNILCat + 11),':',  char(double('B')+ hor -1), num2str(13 + numILCat + numCIFCat + numCOFCat + 2 * numNILCat + 10)));
   [runOffICashIF, ~,  ~] = xlsread(liabilityScenarioFileName, 'Run Off', strcat('B13:CC', num2str(12 + numCIFCat)));
   [runOffICashOF, ~,  ~] = xlsread(liabilityScenarioFileName, 'Run Off', strcat('B', num2str(15 + numCIFCat), ':CC', num2str(14 + numCIFCat + numCOFCat)));
   [quantDK, ~,  ~] = xlsread(liabilityScenarioFileName, 'Quantiles', strcat('A8:', char(double('A')+ hor),  num2str(8 +  numQuant - 1)));
   [quantCF, ~,  ~] = xlsread(liabilityScenarioFileName, 'Quantiles', strcat('A', num2str(11 + numQuant),':', char(double('A')+ hor),  num2str(10 + 2 * numQuant)));
   [sanMass, ~,  ~] = xlsread(liabilityScenarioFileName, 'Liabilities', strcat('B',num2str(26 + numILCat + numCIFCat + numCOFCat + 2 * numNILCat), ':', char(double('B')+ hor -1), num2str(26 + numILCat + numCIFCat + numCOFCat + 2 * numNILCat)));
   [costs, ~,  ~] = xlsread(liabilityScenarioFileName, 'Returns', strcat('D',num2str(12 + 4 * numPsLiab), ':', char(double('D')+ hor -1), num2str(12 + 4 * numPsLiab)));
   
   [spotGov, ~, ~] = xlsread(liabilityScenarioFileName, 'EidgenossenTerminStruktur', 'A4:B14');
   %long end extrapolation
   spotGov = [spotGov;runOff, spotGov(end, 2)];
   discountFactorsDat = [spotGov(:,1),1./((1+spotGov(:,2)).^spotGov(:,1))];
   continousRatesDat = [discountFactorsDat(:,1), -1./discountFactorsDat(:,1).*log(discountFactorsDat(:,2))];
   discountFactorFn = @(t) exp(-t.*(interp1(continousRatesDat(:,1), continousRatesDat(:,2), t, 'linear')));
   
   [~, testList1, ~] = xlsread(liabilityScenarioFileName, 'Run Off', strcat('A13:A',num2str(12 + numTest1)));
   if not(isequal(testList1, iCashIFCat))
      iCashIFCat
      testList1
      msgID = 'createPensionFundLiabilities:badargException';
      msg = 'Incompatibility between insurance cash in flow list in liabilities and run off sheets.';
      throw(MException(msgID,msg));
   end
   
   [~, testList2, ~] = xlsread(liabilityScenarioFileName, 'Run Off', strcat('A', num2str(15 + numTest1), ':A', num2str(14 + numTest1 + numTest2)));
   if not(isequal(testList2, iCashOFCat))
      iCashOFCat
      testList2
      msgID = 'createPensionFundLiabilities:badargException';
      msg = 'Incompatibility between insurance cash out flow list in liabilities and run off sheets.';
      throw(MException(msgID,msg));
   end
   
   [~, testList3, ~] = xlsread(liabilityScenarioFileName, 'Returns', strcat('B10:B', num2str(9 + 4 * numTest3)));
   testList3 = testList3(not(cellfun(@(a)strcmp(a,''), testList3)));
   testList3 = testList3(numTest3 - numNILCat +1:end);
   if not(isequal(testList3, nILCat))
      nILCat
      testList3
      msgID = 'createPensionFundLiabilities:badargException';
      msg = 'Incompatibility between non insurance liability list in liabilities and return sheets.';
      throw(MException(msgID,msg));
   end

result = PensionFundLiabilities(iLCat, nILCat, nILCfRelCat, ilDat, nilDat, iCashIFCat, iCashOFCat, ...
                                iCashIFDat, iCashOFDat, psCat, psDat, vfRes, mvDev, runOffICashIF,...
                                runOffICashOF, quantDK, quantCF, sanMass, costs, discountFactorFn);

end







