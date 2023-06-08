classdef PensionFundLiabilities
    %PensionFundAssets represents the liability side of a pension fund
    
    properties
        
      %Insurance liabilities categories (a cell array of strings)  
      insuranceLiabCat;
      
      %Non insurance liabilities categories (a cell array of strings)
      nonInsuranceLiabCat;
      
      %Cash flow non insurance liabilities categories (a cell array of strings)
      cashFlowRelevantNonInsuranceLiabCat;
      
      %Insurance liabilities data (a numerical vertical vector)
      iLiabDat;
      
      %Non insurance liabilities data (a numerical matrix with two columns)
      niLiabDat;
      
      %Cash inflow insurance liabilities categories (a cell array of strings)
      iCashInFlowsCat;
      
      %Cash outflow insurance liabilities categories (a cell array of strings)
      iCashOutFlowsCat;
      
      %Cash inflow insurance liabilities data (a numerical matrix)
      iCashInFlowsDat;
      
      %Cash outflow insurance liabilities data (a numerical matrix)
      iCashOutFlowsDat;
      
      %Pseudo liabilities categories (a cell array of strings)
      pseudoLiabCat;
      
      %Pseudo liabilities data (a numerical matrix)
      pseudoLiabDat;
	  
      %Initial value fluctuation reserve (a horizontal 2 numerical vector)
      valFluctRes;
      
      %Future market value evolution of non insurance liabilities (a numerical matrix)
      mvDevel;
      
      %Run off cash inflows insurance liabilities (a numerical matrix)
      runOffICashInFlows;
      
      %Run off cash outflows insurance liabilities (a numerical matrix)
      runOffICashOutFlows;
      
      %Quantiles of insurance reserves (a numerical matrix)
      quantileReserves;
      
      %Quantiles of insurance net cash flows (a numerical matrix)
      quantileICashFlows; 
      
      %Non insurance liability costs (a horizontal numerical vector)
      liabilityHardCfCosts;
      
      %Pension fund administration costs (a horizontal numerical vector)
      administrationCosts;
      
      %Discount factor function given as government term structure (a function)
      discountFactorFn
      
      %Time buckets for insurance liability cash flows replication as gov
      %bonds portfolio
      partitionGovBondMat = [0, 1; 2, 5; 6, 10; 11, 15; 16, inf];
      
    end
    
    %public methods
    methods (Access = public)
        
      %Constructor
      function  obj = PensionFundLiabilities(iLCat0, nILCat0, nILCfRelCat0, ilDat0, nilDat0, iCashIFCat0,... 
	                                         iCashOFCat0, iCashIFDat0, iCashOFDat0, psCat0, psDat0, vfl0,...
	                                         mvDev0, runOffICashIF0, runOffICashOF0, qRes0, qCF0, hcfc0, admCosts0, discountFactorFn0)
       
       obj.insuranceLiabCat = iLCat0;
       obj.nonInsuranceLiabCat = nILCat0;
       obj.cashFlowRelevantNonInsuranceLiabCat = nILCfRelCat0;
       obj.iLiabDat = ilDat0;
       obj.niLiabDat = nilDat0;
       obj.iCashInFlowsCat = iCashIFCat0;
       obj.iCashOutFlowsCat = iCashOFCat0;
       obj.iCashInFlowsDat = iCashIFDat0;
       obj.iCashOutFlowsDat = iCashOFDat0;
       obj.pseudoLiabCat = psCat0;
       obj.pseudoLiabDat = psDat0;
       obj.valFluctRes = vfl0;
       obj.mvDevel = mvDev0;
       obj.runOffICashInFlows = runOffICashIF0;
       obj.runOffICashOutFlows = runOffICashOF0;
      
       obj.quantileReserves = qRes0;
       obj.quantileICashFlows = qCF0;
       
       obj.liabilityHardCfCosts = hcfc0;%Sanierungsmassnahmen
       obj.administrationCosts = admCosts0;
       obj.discountFactorFn = discountFactorFn0;
   
      end
      
      function  iLCat0 = getInsuranceLiabilitiesCategories(obj)
         iLCat0 = obj.insuranceLiabCat;            
      end
      
      function  discountFactorFn0 = getGovDiscountFactorFunction(obj)
         discountFactorFn0 = obj.discountFactorFn;            
      end
      
      function  mvfl0 = getInitialValueFluctuationReserveValue(obj)
         vfl0 = obj.valFluctRes;   
         mvfl0 = vfl0(1);
      end
      
      function  hcfc0 = getLHardCfCosts(obj)
         %Sanierungsmassnahmen and other liability pension fund costs not covered by insurance liability administration costs 
         %e.g. administration costs non insurance liabilities
         hcfc0 = obj.liabilityHardCfCosts;            
      end
      
      function  admCosts0 = getLAdministrationCosts(obj)
         %Pension fund insurance liability administration costs as hard cash flows
         admCosts0 = obj.administrationCosts;            
      end
      
      function  qRes0 = getQuantilesReserves(obj)
         qRes0 = obj.quantileReserves;            
      end
      
      function  qCF0 = getQuantilesInsuranceCashFlows(obj)
         qCF0 = obj.quantileICashFlows;            
      end
      
      function  nILCat0 = getNonInsuranceLiabilitiesCategories(obj)
         nILCat0 = obj.nonInsuranceLiabCat;          
      end
      
      function  iCashIFCat0 = getInsuranceLiabilitiesCashInFlowsCategories(obj)
          iCashIFCat0 = obj.iCashInFlowsCat;           
      end
      
      function  iCashOFCat0 = getInsuranceLiabilitiesCashOutFlowsCategories(obj)
         iCashOFCat0 = obj.iCashOutFlowsCat;            
      end
      
      function  psCat0 = getPseudoLiabilitiesCategories(obj)
          psCat0 = obj.pseudoLiabCat;           
      end
     
      function  result = isInsuranceCashInFlowPensionFundLiabilityQ(obj, cfliab)
          result = ismember(cfliab, obj.iCashInFlowsCat);           
      end
      
      function  nILCfRelCat0 = getNonInsuranceCashFlowRelevantLiabilitiesCategories(obj)
          nILCfRelCat0 = obj.cashFlowRelevantNonInsuranceLiabCat;            
      end
      
      function  result = isInsurancePensionFundLiabilityQ(obj, liab0)
         result = ismember(liab0, obj.insuranceLiabCat);           
      end
      
      function  result = isNonInsurancePensionFundLiabilityQ(obj, liab0)
         result = ismember(liab0, obj.nonInsuranceLiabCat);           
      end
      
      function result = isInsuranceCashInflowPensionFundLiabilityQ(obj, liab0)
         result = ismember(liab0, obj.iCashInFlowsCat);           
      end
      
      function  result = isInsuranceCashOutFlowPensionFundLiabilityQ(obj, liab0)
         result = ismember(liab0, obj.iCashOutFlowsCat);           
      end
      
      function  result = isPseudoPensionFundLiabilityQ(obj, liab0)
         result = ismember(liab0, obj.pseudoLiabCat);           
      end
      
      function  result = getInsuranceLiabilityValues(obj, liab0)
        if isInsurancePensionFundLiabilityQ(obj, liab0)
         iLiabDat0 = obj.iLiabDat;
         result = iLiabDat0(insuranceLiabilityNum(obj, liab0), :);
        else
         msgID = 'getInsuranceLiabilityValues:badargException';
         msg = strcat(liab0, ' not in insurance liabilities category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getNonInsuranceLiabilityInitialValue(obj, liab0)
        if isNonInsurancePensionFundLiabilityQ(obj, liab0)
         niLiabDat0 = obj.niLiabDat;
         result = niLiabDat0(nInsuranceLiabilityNum(obj, liab0), 1);
        else
         msgID = 'getNonInsuranceLiabilityInitialValue:badargException';
         msg = strcat(liab0, ' not in non insurance liabilities category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      
      function  result = getNonInsuranceLiabilityValuePrescribedEvolution(obj, liab0)
        if isNonInsurancePensionFundLiabilityQ(obj, liab0)
         mvDevel0 = obj.mvDevel;
         result = mvDevel0(nInsuranceLiabilityNum(obj, liab0), :);
        else
         msgID = 'getNonInsuranceLiabilityValuePrescribedEvolution:badargException';
         msg = strcat(liab0, ' not in non insurance liabilities category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getNonInsuranceLiabilityWeightPrescribedEvolution(obj, liab0)
        if isNonInsurancePensionFundLiabilityQ(obj, liab0)
         mvDevel0 = obj.mvDevel;
         result = mvDevel0(nInsuranceLiabilityNum(obj, liab0), :) ./ sum(mvDevel0);
        else
         msgID = 'getNonInsuranceLiabilityWeightPrescribedEvolution:badargException';
         msg = strcat(liab0, ' not in non insurance liabilities category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getNonInsuranceLiabilityInitialWeight(obj, liab0)
        if isNonInsurancePensionFundLiabilityQ(obj, liab0)
         niLiabDat0 = obj.niLiabDat;
         result = niLiabDat0(nInsuranceLiabilityNum(obj, liab0), 1) ./ sum(niLiabDat0(:,1));
        else
         msgID = 'getNonInsuranceLiabilityInitialWeight:badargException';
         msg = strcat(liab0, ' not in non insurance liabilities category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result =  getTotalInsuranceRunOffCashFlows(obj)
         inflows = sum(cellfunVectorOutput(@ (arg) getInsuranceLiabilityRunOffCashInFlows(obj,arg), getInsuranceLiabilitiesCashInFlowsCategories(obj)), 1);
         outflows = sum(cellfunVectorOutput(@ (arg) getInsuranceLiabilityRunOffCashOutFlows(obj,arg), getInsuranceLiabilitiesCashOutFlowsCategories(obj)), 1);
         result = outflows - inflows;
      end
      
      function  result =  getTotalInsuranceCashFlows(obj)
         inflows = sum(cellfunVectorOutput(@ (arg) getInsuranceLiabilityCashInFlows(obj,arg), getInsuranceLiabilitiesCashInFlowsCategories(obj)), 1);
         outflows = sum(cellfunVectorOutput(@ (arg) getInsuranceLiabilityCashOutFlows(obj,arg), getInsuranceLiabilitiesCashOutFlowsCategories(obj)), 1);
         result = outflows - inflows;
      end
      
      function result = getInsuranceLiabilityInitialMarketWeights(obj)
        govCat = setdiff(getPseudoLiabilitiesCategories(obj), getNonInsuranceLiabilitiesCategories(obj));
        cfs = getTotalInsuranceRunOffCashFlows(obj);
        den = getTotalInsuranceLiabilityInitialMarketValue(obj);
        timeBuckets = obj.partitionGovBondMat;
        disc =  getGovDiscountFactorFunction(obj);
        times = 1:length(cfs);
        ds = disc(times);
        nGovCat = length(govCat);
        pvs = zeros(nGovCat,1); 
        for govCatNum = 1:nGovCat 
            pos = times(and(times>=timeBuckets(govCatNum, 1), times<=timeBuckets(govCatNum, 2)));
            pvs(govCatNum) = sum(cfs(pos).* ds(pos));
        end
        result = pvs/den;
      end    
      
      function  result = getInsuranceLiabilityCashInFlows(obj, liab0)
        if isInsuranceCashInFlowPensionFundLiabilityQ(obj, liab0)
         iCashInFlowsDat0 = obj.iCashInFlowsDat;
         result = iCashInFlowsDat0(insuranceLiabilityCashInFlowNum(obj, liab0), :);
        else
         msgID = 'getInsuranceLiabilityCashInFlows:badargException';
         msg = strcat(liab0, ' not in insurance cash inflows category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getInsuranceLiabilityCashOutFlows(obj, liab0)
        if isInsuranceCashOutFlowPensionFundLiabilityQ(obj, liab0)
         iCashOutFlowsDat0 = obj.iCashOutFlowsDat;
         result = iCashOutFlowsDat0(insuranceLiabilityCashOutFlowNum(obj, liab0), :);
        else
         msgID = 'getInsuranceLiabilityCashOutFlows:badargException';
         msg = strcat(liab0, ' not in insurance cash outflows category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getInsuranceLiabilityRunOffCashInFlows(obj, liab0)
        if isInsuranceCashInFlowPensionFundLiabilityQ(obj, liab0)
         runOffICashInFlows0 = obj.runOffICashInFlows;
         result = runOffICashInFlows0(insuranceLiabilityCashInFlowNum(obj, liab0), :);
        else
         msgID = 'getInsuranceLiabilityRunOffCashInFlows:badargException';
         msg = strcat(liab0, ' not in insurance cash inflows category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function  result = getInsuranceLiabilityRunOffCashOutFlows(obj, liab0)
        if isInsuranceCashOutFlowPensionFundLiabilityQ(obj, liab0)
         runOffICashOutFlows0 = obj.runOffICashOutFlows;
         result = runOffICashOutFlows0(insuranceLiabilityCashOutFlowNum(obj, liab0), :);
        else
         msgID = 'getInsuranceLiabilityRunOffCashOutFlows:badargException';
         msg = strcat(liab0, ' not in insurance cash outflows category list.');
         throw(MException(msgID,msg)); 
        end              
      end
      
      function result = getTotalInsuranceLiabilityInitialMarketValue(obj)
          iRunOffCat0 = getInsuranceLiabilitiesCashInFlowsCategories(obj);
          oRunOffCat0 = getInsuranceLiabilitiesCashOutFlowsCategories(obj);
          cfs = sum(cellfunVectorOutput(@(arg) getInsuranceLiabilityRunOffCashOutFlows(obj, arg), oRunOffCat0), 1)...
                - sum(cellfunVectorOutput(@(arg) getInsuranceLiabilityRunOffCashInFlows(obj, arg), iRunOffCat0), 1);
          disc =  getGovDiscountFactorFunction(obj);
          times = 1:length(cfs);
          result = sum(cfs.*disc(times));
      end
      
      function  result = getPseudoLiabilityExpectedTotalReturns(obj, liab0)
        if isPseudoPensionFundLiabilityQ(obj, liab0)
           pseudoLiabDat0 = obj.pseudoLiabDat;
           result = pseudoLiabDat0(4 * (pseudoLiabNum(obj, liab0)-1) + 3, :);
         else
           msgID = 'getPseudoLiabilityExpectedTotalReturns:badargException';
           msg = strcat(liab0, ' not in pseudo liabilities category list.');
           throw(MException(msgID,msg)); 
         end              
      end
      
      function  result = getCumulativePseudoLiabilityExpectedTotalReturns(obj, liab0)
        if isPseudoPensionFundLiabilityQ(obj, liab0)
           rets = getPseudoLiabilityExpectedTotalReturns(obj, liab0);
           result = cumprod(1+rets)-1; 
        else
           msgID = 'getCumulativePseudoLiabilityExpectedTotalReturns:badargException';
           msg = strcat(liab0, ' not in pseudo liabilities  category list.');
           throw(MException(msgID,msg)); 
        end              
      end
     
      function  result = getPseudoLiabilityExpectedPriceReturns(obj, liab0)
       if isPseudoPensionFundLiabilityQ(obj, liab0)
          pseudoLiabDat0 = obj.pseudoLiabDat;
          result = pseudoLiabDat0(4 * (pseudoLiabNum(obj, liab0)-1) + 1, :); 
       else
          msgID = 'getPseudoLiabilityExpectedPriceReturns:badargException';
          msg = strcat(liab0, ' not in pseudo liabilities category list.');
          throw(MException(msgID,msg)); 
       end                  
      end
     
      function  result = getPseudoLiabilityExpectedRunningYields(obj, liab0)
        if isPseudoPensionFundLiabilityQ(obj, liab0)
           pseudoLiabDat0 = obj.pseudoLiabDat;
           result = pseudoLiabDat0(4 * (pseudoLiabNum(obj, liab0)-1) + 2, :); 
        else
           msgID = 'getPseudoLiabilityExpectedRunningYields:badargException';
           msg = strcat(liab0, ' not in pseudo liabilities category list.');
           throw(MException(msgID,msg)); 
        end                     
      end
     
      function  result = getPseudoLiabilityTotalReturnVolatilities(obj, liab0)
       if isPseudoPensionFundLiabilityQ(obj, liab0)
          pseudoLiabDat0 = obj.pseudoLiabDat;
          result = pseudoLiabDat0(4 * (pseudoLiabNum(obj, liab0)-1) + 4, :); 
       else
         msgID = 'getPseudoLiabilityTotalReturnVolatilities:badargException';
         msg = strcat(liab0, ' not in pseudo liabilities category list.');
         throw(MException(msgID,msg)); 
       end                     
      end
     
      function  result = getCumulativePseudoLiabilityTotalReturnVolatilities(obj, liab0)
       if isPseudoPensionFundLiabilityQ(obj, liab0)
          rets = getCumulativePseudoLiabilityExpectedTotalReturns(obj, liab0);
          vols = getPseudoLiabilityTotalReturnVolatilities(obj, liab0);
          cumRets = getCumulativePseudoLiabilityExpectedTotalReturns(obj, liab0); 
          T = length(vols);
          vars = zeros(1,T);
          vars(1) = (vols(1))^2;
          for t=2:T
            vars(t) =  vars(t-1) * (vols(t)^2 + (1 + rets(t))^2) +  (1 + cumRets(t-1))^2 * vols(t)^2; 
          end
          result = sqrt(vars);
       else
         msgID = 'getCumulativePseudoLiabilityTotalReturnVolatilities:badargException';
         msg = strcat(liab0, ' not in  pseudo liabilities category list.');
         throw(MException(msgID,msg)); 
       end                          
      end
      
      function  cf = getExpectedNonInsuranceLiabilityCashFlows(obj, liab0)
       if isNonInsurancePensionFundLiabilityQ(obj, liab0)
          rets = getPseudoLiabilityExpectedTotalReturns(obj, liab0);
          T = length(rets);
          pn = zeros(1, T + 1);
          cf = zeros(size(pn));
          presLevel = getNonInsuranceLiabilityValuePrescribedEvolution(obj, liab0);
          pn(1) = getNonInsuranceLiabilityInitialValue(obj, liab0);
          pn(2:end) = presLevel;
          for t=2:T+1
            cf(t) = pn(t - 1) * (1 + rets(t - 1)) - presLevel(t - 1); 
          end
          cf(not(isnumeric(cf))) = 0;
          
       else
         msgID = 'getExpectedNonInsuranceLiabilityCashFlows:badargException';
         msg = strcat(liab0, ' not in non isurance liabilities category list.');
         throw(MException(msgID,msg)); 
       end                          
      end
     
    end  
    
    %private methods
    methods (Access = private)
       
     function  result = pseudoLiabNum(obj, liab0)
         result = find(ismember(obj.pseudoLiabCat, liab0));            
     end
     
     function  result = insuranceLiabilityNum(obj, liab0)
         result = find(ismember(obj.insuranceLiabCat, liab0));            
     end
     
     function  result = nInsuranceLiabilityNum(obj, liab0)
         result = find(ismember(obj.nonInsuranceLiabCat, liab0));            
     end
     
     function  result = insuranceLiabilityCashInFlowNum(obj, liab0)
         result = find(ismember(obj.iCashInFlowsCat, liab0));            
     end
    
     function  result = insuranceLiabilityCashOutFlowNum(obj, liab0)
         result = find(ismember(obj.iCashOutFlowsCat, liab0));            
     end
     
    end  
    
end

   

