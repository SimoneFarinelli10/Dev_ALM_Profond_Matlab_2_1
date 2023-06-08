function optStruct = prepareOptionsCalcPensionFund(outputFolderName)
  
  optStruct = struct();
   
  %Asset cash flows inforce strategy   
  optStruct.('inForceAssCfCsvFileName') = [outputFolderName, 'AssetCashFlowsInForce.csv'];
  
  %Asset cash flows alternative strategy 
  optStruct.('altAssCfCsvFileName') = [outputFolderName, 'AssetCashFlowsAlternative.csv'];
        
  %Insurance cash flows
  optStruct.('iCfCsvFileName') = [outputFolderName, 'InsuranceCashFlows.csv']; 
  
  %Non Insurance cash flows
  optStruct.('nICfCsvFileName') = [outputFolderName, 'NonInsuranceCashFlows.csv'];
  
  %Quantiles reserve
  optStruct.('quantDKPlotJpgFileName') = [outputFolderName, 'QuantilesDK.jpg']; 
  optStruct.('DKCsvFileName') = [outputFolderName, 'QuantilesDK.csv']; 
  
  %Quantiles insurance cashflows
  optStruct.('quantCFPlotJpgFileName') = [outputFolderName, 'QuantilesInsuranceCF.jpg'];
  optStruct.('CFCsvFileName') = [outputFolderName, 'QuantilesInsuranceCF.csv'];
  
  %Cumulated total returns for inforce asset strategy
  optStruct.('inForceCumTotalRetDistJpgFileName') = [outputFolderName, 'DistributionsCumulatedTotalReturnsInForce.jpg']; 
  optStruct.('inForceCumTotalRetQuantJpgFileName') = [outputFolderName, 'QuantilesCumulatedTotalReturnsInForce.jpg']; 
  optStruct.('inForceCumTotalRetCsvFileName') = [outputFolderName, 'QuantilesCumulatedTotalReturnsInForce.csv'];
  
  %Total returns for inforce asset strategy
  optStruct.('inForceTotalRetDistJpgFileName') = [outputFolderName, 'DistributionsTotalReturnsInForce.jpg']; 
  optStruct.('inForceTotalRetQuantJpgFileName') = [outputFolderName, 'QuantilesTotalReturnsInForce.jpg']; 
  optStruct.('inForceTotalRetCsvFileName') = [outputFolderName, 'QuantilesTotalReturnsInForce.csv']; 
  
  %Reported funding ratio for inforce asset strategy
  optStruct.('inForceRfrDistJpgFileName') = [outputFolderName, 'DistributionsRegulatoryFundingRatioInForce.jpg']; 
  optStruct.('inForceRfrQuantJpgFileName') = [outputFolderName, 'QuantilesRegulatoryFundingRatioInForce.jpg']; 
  optStruct.('inForceRfrCsvFileName') = [outputFolderName, 'QuantilesRegulatoryFundingRatioInForce.csv'];
  
  %Cumulated total returns for alternative asset strategy
  optStruct.('altCumTotalRetDistJpgFileName') = [outputFolderName, 'DistributionsCumulatedTotalReturnsAlternative.jpg'];
  optStruct.('altCumTotalRetQuantJpgFileName') = [outputFolderName, 'QuantilesCumulatedTotalreturnsAlternative.jpg']; 
  optStruct.('altCumTotalRetCsvFileName') = [outputFolderName, 'QuantilesCumulatedTotalReturnsAlternative.csv']; 
  
  %Total returns for alternative asset strategy
  optStruct.('altTotalRetDistJpgFileName') = [outputFolderName, 'DistributionsTotalReturnsAlternative.jpg'];
  optStruct.('altTotalRetQuantJpgFileName') = [outputFolderName, 'QuantilesTotalReturnsAlternative.jpg']; 
  optStruct.('altTotalRetCsvFileName') = [outputFolderName, 'QuantilesTotalReturnsAlternative.csv']; 
  
  %Reported funding ratio for alternative asset strategy
  optStruct.('altRfrDistJpgFileName') = [outputFolderName, 'DistributionsRegulatoryFundingRatioAlternative.jpg'];
  optStruct.('altRfrQuantJpgFileName') = [outputFolderName, 'QuantilesRegulatoryFundingRatioAlternative.jpg']; 
  optStruct.('altRfrCsvFileName') = [outputFolderName, 'QuantilesRegulatoryFundingRatioAlternative.csv']; 
  
  %Shortfall probabilities for inforce asset strategy
  optStruct.('inForceProbJpgFileName') = [outputFolderName, 'ShortfallProbabilitiesInForce.jpg'];
  optStruct.('inForceProbCsvFileName') = [outputFolderName, 'ShortfallProbabilitiesInForce.csv'];
  
  %Expected Shortfall Funding Ratio for inforce asset strategy
  optStruct.('inForceESJpgFileName') = [outputFolderName, 'ESInForce.jpg'];
  optStruct.('inForceESCsvFileName') = [outputFolderName, 'ESInForce.csv'];
  
  %Shortfall probabilities for alternative asset strategy
  optStruct.('altProbJpgFileName') = [outputFolderName, 'ShortfallProbabilitiesAlternative.jpg'];
  optStruct.('altProbCsvFileName') = [outputFolderName, 'ShortfallProbabilitiesAlternative.csv'];
  
  %Expected Shortfall Funding Ratio for alternative asset strategy
  optStruct.('altESJpgFileName') = [outputFolderName, 'ESAlternative.jpg'];
  optStruct.('altESCsvFileName') = [outputFolderName, 'ESAlternative.csv'];
  
  %Initial regulatory funding ratio
  optStruct.('initialReportedFundingRatioCsvFileName') = [outputFolderName, 'InitialRegulatoryFundingRatio.csv']; 
  
  %ALM cumulated total return efficient frontier
  optStruct.('optPortCsvFileName') = [outputFolderName, 'ALMEfficientPortfolios.csv'];
  optStruct.('effFrontJpgFileName') = [outputFolderName, 'ALMEfficientFrontier.jpg'];

  %Target value fluctuation reserves
  optStruct.('inForceTargetReserveCsvFileName') = [outputFolderName, 'TargetReserveInForce.csv'];
  optStruct.('inForceCumulatedTargetReserveCsvFileName') = [outputFolderName, 'CumulatedTargetReserveInForce.csv'];
  optStruct.('altTargetReserveCsvFileName') = [outputFolderName, 'TargetReserveAlternative.csv']; 
  optStruct.('altCumulatedTargetReserveCsvFileName') = [outputFolderName, 'CumulatedTargetReserveAlternative.csv']; 

  %Asset strategies
  optStruct.('inForceAlloAssetJpgFileName') = [outputFolderName, 'AlloAssetInForce.jpg']; 
  optStruct.('inForceAlloAssetCsvFileName') = [outputFolderName, 'AlloAssetInForce.csv']; 
  optStruct.('altAlloAssetJpgFileName') = [outputFolderName, 'AlloAssetAlternative.jpg']; 
  optStruct.('altAlloAssetCsvFileName') = [outputFolderName, 'AlloAssetAlternative.csv']; 
  
  %Balance sheet
  optStruct.('inForceBsCsvFileName') = [outputFolderName, 'BsInForce.csv']; 
  optStruct.('inForceBsJpgFileName') = [outputFolderName, 'BsInForce.jpg']; 
  optStruct.('altBsCsvFileName') = [outputFolderName, 'BsAlternative.csv']; 
  optStruct.('altBsJpgFileName') = [outputFolderName, 'BsAlternative.jpg']; 
  
         
  
end



