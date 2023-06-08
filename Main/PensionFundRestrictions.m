classdef PensionFundRestrictions
     %PensionFundRestrictions represents the restrictions on the asset side of a pension fund
    
    properties
       
       %Asset categories (a vertical cell array of strings) 
       assetCat;
       
       %Restriction coefficients (a numerical matrix)
       matrixDat;
       
       
    end
    
    %public methods
    methods (Access = public)
        
     %Constructor
     function  obj = PensionFundRestrictions(assetCat0, matrixDat0)
      obj.assetCat = assetCat0;
      obj.matrixDat = matrixDat0;         
     end  
        
     %Get methods
     function  assetCat0 = getAssetsCategories(obj)
      assetCat0 = obj.assetCat;            
     end
     
     function  result = getLowerBound(obj, asset0)
       if isAssetRestrictionsPensionFundQ(obj, asset0)
          matDat0 = obj.matrixDat;
          result = matDat0(assetNum(obj, asset0), 1);
       else
          msgID = 'getLowerBound:badargException';
          msg = strcat(asset0, ' not in asset category list.');
          throw(MException(msgID,msg)); 
      end                 
     end
     
     function  result = getUpperBound(obj, asset0)
       if isAssetRestrictionsPensionFundQ(obj, asset0)
          matDat0 = obj.matrixDat;
          result = matDat0(assetNum(obj, asset0), 2);
       else
          msgID = 'getUpperBound:badargException';
          msg = strcat(asset0, ' not in asset category list.');
          throw(MException(msgID,msg)); 
      end                 
     end
     
     function  result = getNonElementaryConstraints(obj)
      matDat0 = obj.matrixDat; 
      result = matDat0(:, 3:end);
     end
     
     function  result = isAssetRestrictionsPensionFundQ(obj, asset0)
      result = ismember(asset0, obj.assetCat);            
     end
     
    end
    
    %private methods
    methods (Access = private)
        
      function  result = assetNum(obj, asset0)
        result = find(ismember(obj.assetCat, asset0));            
      end 
     
    end
    
end




