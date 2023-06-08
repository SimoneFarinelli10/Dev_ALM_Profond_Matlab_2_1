classdef PensionFundAssets < handle
    %PensionFundAssets represents the asset side of a pension fund
    
    properties
      
      %Asset categories (a vertical cell array of strings)
      assetCat;
      
      %Asset data: expected price returns, expected running yields,
      %expected total returns, total return volatilities (a numerical matrix)
      assetData; 
      
      %Asset initial market value (a numerical vertical vector)
      marketV0; 
      
      %Asset future market values (a numerical matrix)
      marketValueEvolution; 
      
      %Asset initial reported value (a numerical vertical vector)
      %reportedV0;
      
      %Alternative Strategies: for every alternative strategy initial weights and asset future market
      %values (a structure, field names are strategy names, fields are
      %numerical matrices, NaN implies relative rebalancing) 
      alternativeStrategies; 
      
      %Securities (a vertical cell array of strings)
      securitiesCat;
      
      %Securities administration costs as percentage of the security wealth
      %at year's beginning (a number --> a vertical numerical vector)
      administrationCosts; 
      
      %Hard cash flows for generic costs, typically for real estate (a
      %horizontal numerical vector)
      assetHardCfCosts;  
      
    end
    
    %public methods
    methods (Access = public)
      
     %Constructor
     function  obj = PensionFundAssets(assetCat0, assetData0, marketV00, marketValueEvolution0,... 
                                       alternativeStrategies0, securitiesCat0, administrationCosts0, assetHardCfCosts0)
       
      obj.assetCat = assetCat0;
      obj.assetData = assetData0;
      obj.marketV0 = marketV00;  
      obj.marketValueEvolution = marketValueEvolution0;  
      obj.alternativeStrategies = alternativeStrategies0;  
      obj.securitiesCat = securitiesCat0;  
      obj.administrationCosts = administrationCosts0;  
      obj.assetHardCfCosts = assetHardCfCosts0;  
               
     end  
        
     %Get/Set methods
     function  assetCat0 = getAssetsCategories(obj)
      assetCat0 = obj.assetCat;            
     end
     
     function  securitiesCat0 = getSecuritiesCategories(obj)
      securitiesCat0 = obj.securitiesCat;            
     end
      
     function  result = getAssetsWithPrescribedFutureValuesCategories(obj, t)
      mat0 = obj.marketValueEvolution;
      hor0  = size(mat0, 2);
      if t <= hor0
       assetCat0 = getAssetsCategories(obj);
       result = assetCat0(not(isnan(mat0(:,t)))); 
      else
       msgID = 'getAssetsWithPrescribedFutureValuesCategories:badargException';
       msg = ['Time ', num2str(t), ' is beyond final horizon ', num2str(hor0), '.'];
       throw(MException(msgID,msg)); 
      end
     end
     
     function  result = getAssetPrescribedFutureValues(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
         mve = obj.marketValueEvolution;
         result = mve(assetNum(obj, asset0), :);
      else
         msgID = 'getAssetPrescribedFutureValues:badargException';
         msg = strcat(asset0, ' not in asset category list.');
         throw(MException(msgID,msg)); 
      end    
     end
    
     function  result = getAssetsWithoutPrescribedFutureValuesCategories(obj, t)
      mat0 = obj.marketValueEvolution;
      hor0  = size(mat0, 2);
      if t <= hor0
       assetCat0 = getAssetsCategories(obj);
       result = assetCat0(isnan(mat0(:,t))); 
      else
       msgID = 'getAssetsWithoutPrescribedFutureValuesCategories:badargException';
       msg = ['Time ', num2str(t), ' is beyond final horizon ', num2str(hor0), '.'];
       throw(MException(msgID,msg)); 
      end           
     end
     
     function  result = getAAdministrationCosts(obj)
      result = obj.administrationCosts;            
     end
     
     function  result = getAHardCfCosts(obj)
      result = obj.assetHardCfCosts;            
     end
     
     function  result = getAssetExpectedTotalReturns(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       mat0 = obj.assetData;
       result = mat0(4 * (assetNum(obj, asset0)-1) + 3, :); 
      else
       msgID = 'getAssetExpectedTotalReturns:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end              
     end
     
     function  result = getCumulativeAssetExpectedTotalReturns(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       rets = getAssetExpectedTotalReturns(obj, asset0);
       result = cumprod(1+rets)-1; 
      else
       msgID = 'getCumulativeAssetExpectedTotalReturns:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end              
     end
     
     function  result = getAssetExpectedPriceReturns(obj, asset0)
       if isAssetPensionFundAssetQ(obj, asset0)
       mat0 = obj.assetData;
       result = mat0(4 * (assetNum(obj, asset0)-1) + 1, :); 
      else
       msgID = 'getAssetExpectedPriceReturns:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                  
     end
     
     function  result = getAssetExpectedRunningYields(obj, asset0)
       if isAssetPensionFundAssetQ(obj, asset0)
       mat0 = obj.assetData;
       result = mat0(4 * (assetNum(obj, asset0)-1) + 2, :); 
      else
       msgID = 'getAssetExpectedRunningYields:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                     
     end
     
     function  result = getAssetTotalReturnVolatilities(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       mat0 = obj.assetData;
       result = mat0(4 * (assetNum(obj, asset0)-1) + 4, :); 
      else
       msgID = 'getAssetTotalReturnVolatilities:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                 
     end
     
     function  result = getCumulativeAssetTotalReturnVolatilities(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       rets = getAssetExpectedTotalReturns(obj, asset0);
       vols = getAssetTotalReturnVolatilities(obj, asset0);
       cumRets = getCumulativeAssetExpectedTotalReturns(obj, asset0); 
       T = length(vols);
       vars = zeros(1,T);
       vars(1) = (vols(1))^2;
       for t=2:T
         vars(t) =  vars(t-1) * (vols(t)^2 + (1 + rets(t))^2) +  (1 + cumRets(t-1))^2 * vols(t)^2; 
       end
       result = sqrt(vars);
      else
       msgID = 'getCumulativeAssetTotalReturnVolatilities:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                          
     end
      
     function  result = isAssetPensionFundAssetQ(obj, asset0)
      result = ismember(asset0, obj.assetCat);            
     end
          
     function  result = getAssetInitialMarketValue(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       v0 = obj.marketV0;
       result = v0(assetNum(obj, asset0)); 
      else
       msgID = 'getAssetInitialMarketWeight:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                      
     end
     
     function  result = getAssetInitialMarketWeight(obj, asset0)
      if isAssetPensionFundAssetQ(obj, asset0)
       v0 = obj.marketV0;
       result = v0(assetNum(obj, asset0))/sum(v0); 
      else
       msgID = 'getAssetInitialMarketWeight:badargException';
       msg = strcat(asset0, ' not in asset category list.');
       throw(MException(msgID,msg)); 
      end                      
     end
      
     function  result = getAlternativeStrategies(obj)
      result = obj.alternativeStrategies;            
     end
     
     function  setAHardCfCosts(obj, costs0)
      hor0 = length(obj.assetHardCfCosts);
      if and(length(costs0) == hor0, isnumeric(costs0))
       obj.assetHardCfCosts = costs0; 
      else
       msgID = 'setAHardCfCosts:badargException';
       msg = ['New cost vector is not acceptable.'];
       throw(MException(msgID,msg)); 
      end                       
     end
    
     function  setAssetPrescribedFutureValues(obj, asset0, values0)
      mat0 = obj.marketValueEvolution;
      hor0  = size(mat0, 2);
      if not(isAssetPensionFundAssetQ(obj, asset0))
        msgID = 'setAssetPrescribedFutureValues:badargException';
        msg = strcat(asset0, ' not in asset category list.');
        throw(MException(msgID,msg)); 
      elseif not(length(values0)==hor0)
        msgID = 'setAssetPrescribedFutureValues:badargException';
        msg = ['Length of new cost vector not compatible with horizon ', num2str(hor0),'.'];
        throw(MException(msgID,msg));   
      else    
       pos = assetNum(obj, asset0);
       mat0(pos,:) = values0;
       obj.marketValueEvolution = mat0;
      end         
     end
         
    end
    
    %private methods
    methods (Access = private)
        
     function  result = assetNum(obj, asset0)  
      result = find(ismember(obj.assetCat, asset0)) ;            
     end  
       
    end    
    
end
