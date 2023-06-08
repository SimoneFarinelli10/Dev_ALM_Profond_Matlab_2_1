classdef PensionFundProjection < handle
    %PensionFundProjection projects the whole balance sheet of a pension
    %fund into the future and allows for a complete ALM study
    
    
    properties
        
        %Asset side (an object of class PensionFundAssets)
        pfass; 
        
        %Liability side (an object of class PensionFundLiabilities)
        pfliab; 
        
        %Asset and liability correlations (a struct)
        corr; 
        
        %Asset restrictions for reporting (an object of class PensionFundRestrictions)
        pfrestr; 
        
        %Asset restrictions for optimization (an object of class PensionFundRestrictions)
        pfrestrOpt; 
        
        %Reporting date (serial date number)
        presentDate;
        
        %Simulated total returns for asset and (pseudo-)liability categories 
        %(a three dimensional numeric array: simulation dimension, time dimension and asset/liability dimension 
        simulTotRet; 
        
        %Simulated pension fund funding ratio (a two dimensional
        %numeric array: simulation dimension and time dimension
        simulFunRat;
          
        %Simulated pension fund asset side (a three dimensional
        %numeric array: simulation dimension,  time dimension and asset dimension 
        assetSide; 
        
        %Simulated pension fund asset side (a three dimensional
        %numeric array: simulation dimension,  time dimension and liability dimension 
        insuranceLiabilitySide; 
        
        %Initial market asset strategy weights (a vertical numeric vector)
        initialAssetStrategy; 
        
        %Market asset strategy weights (a numeric matrix)
        assetStrategy; 
        
        %Asset strategy name (a string)
        assetStrategyName = 'SAA';
        
        %Total administration costs on asset and liability side (a
        %horizontal numeric vector)
        administrationCosts; 
       
        %Optimization inputs: horizon and asset categories to freeze (a struct)
        inputsOpt; 
        
        %numbers of simulation (an integer)
        numSim; 
        
        %numbers of points on the efficient frontier (an integer)
        numOpt
        
    end
    
%     Methods for class PensionFundProjection:
% 
%     PensionFundProjection                              
%     calcAssetMarketStrategy                            
%     catNum                                             
%     checkRestrictionsOverTime                          
%     clearCache                                         
%     getAIncreaseFactor                                 
%     getAReturns                                        
%     getAdministrationCosts                             
%     getAssetCashFlows                                  
%     getAssetStrategy                             
%     getAssetReportedStrategy                           
%     getCategories                                      
%     getCfRelevantPNIncreaseFactor                      
%     getCorrelations                                    
%     getExpectedTotalInsuranceCashFlows                 
%     getExpectedTotalLiabilityCashFlows                 
%     getExpectedTotalNonInsuranceCashFlows              
%     getInitialMarketBasedFundingRatio                  
%     getInitialReportedFundingRatio                     
%     getLHardCashFlowCosts                              
%     getNumSim                                          
%     getNumberOptimalPortfolios                         
%     getTargetValueFluctuationReserves            
%     getCumulatedTargetValueFluctuationReserves                 
%     getPNCashFlows                                     
%     getPNIncreaseFactor                                
%     getPNReturns                                       
%     getPensionFundAssets                               
%     getPensionFundLiabilities                          
%     getPensionFundOptimizationRestrictions             
%     getPensionFundPresentDate                          
%     getPensionFundRestrictions                         
%     getPieChartReportedAssetStrategy                   
%     getRegulatoryBalanceSheetProjection                
%     getReportedFundingRatioProbability                 
%     getRestrictionsFunction                            
%     getSimulatedTotalReturns                           
%     getTotalAssetCashFlows                             
%     getTotalAssetInitialValue                    
%     getTotalAssetInitialReportedValue                  
%     getTotalCashFlows                                  
%     getTotalInsuranceCashFlows                         
%     getTotalInsuranceLiabilityInitialMarketValue       
%     getTotalInsuranceValues                    
%     getTotalInsuranceRunOffCashFlows                   
%     getTotalLiabilityCashFlows                         
%     getTotalNonInsuranceLiabilityInitialMarketValue    
%     getTotalNonInsuranceLiabilityInitialValue  
%     getVCashFlows                                      
%     getVIncreaseFactor                                 
%     getVReturns                                        
%     getaIncreaseFactor                                 
%     initialize                                         
%     isCategoryPensionFundProjectionQ                   
%     projectionCashFlowInsurance                        
%     projectionCumulatedTotalReturnsAssetSide           
%     projectionDefaultProbabilities                     
%     projectionAssetSide                          
%     projectionMarketBasedFundingRatio                  
%     projectionMarketInsuranceLiabilitySide             
%     projectionReportedAssetSide                        
%     projectionFundingRatio                     
%     projectionReportedInsuranceLiabilitySide           
%     projectionTotalReturnsAssetSide                    
%     setNumSim                                          
%     setNumberOptimalPortfolios                         
% 
% 
%     Static methods:
% 
%     barPlot                                            
%     createGenericTableOutput                           
%     createTableOutput                                  
%     getStatistics                                      
%     funRat                                          
%     plotDensity                                        
%     plotPieChartStrategy                               
%     funRat                                          
%     sigma 
%     
    %static methods (local helper functions)
    methods(Static)
    
        function stat = getStatistics(data) 
         stat = [mean(data); std(data, 1); quantile(data,  [0.95, 0.75, 0.50, 0.25, 0.05])];
        end 
        
        function tableOutput = createTableOutput(statValues, text) 
         rowNamesValues = {'Expectation'; 'Volatility'; 'Quantile 0.95'; 'Quantile 0.75'; 'Quantile 0.50'; 'Quantile 0.25'; 'Quantile 0.05'};
         hor = size(statValues,2);
         columnNamesValues = cell(1, hor);
         for t=1:hor
            columnNamesValues(t) = {['Year_', num2str(t)]};
         end    
         tableOutput = array2table(statValues);
         tableOutput.Properties.VariableNames = columnNamesValues;
         tableOutput.Properties.RowNames = rowNamesValues;
         dimNames = tableOutput.Properties.DimensionNames;
         if nargin >= 2
            dimNames{1} = text;
         end
         tableOutput.Properties.DimensionNames = dimNames;
        end
        
        function tableOutput = createGenericTableOutput(numValues, rowNamesValues, columnNamesValues, text) 
         tableOutput = array2table(numValues);
         tableOutput.Properties.VariableNames = columnNamesValues;
         tableOutput.Properties.RowNames = rowNamesValues;
         dimNames = tableOutput.Properties.DimensionNames;
         if nargin >= 4
            dimNames{1} = text;
         end
         tableOutput.Properties.DimensionNames = dimNames;
        end
          
        function barPlot(statValues, titleText, unit, scal, outJpgFileName) 
              
          function result = makeBarChartInput(quantiles)
           result = quantiles - [zeros(1, size(quantiles,2));quantiles(1:end-1,:)];
          end    
          fig = figure;
          plotInput = scal * (makeBarChartInput(statValues(end:-1:3,:)))';
          barChartPlot = bar(plotInput, 'stacked', 'EdgeColor', 'w');
          barChartPlot(1).FaceColor = 'w';
          bl = barChartPlot.BaseLine;
          bl.Visible = 'off';
          plotAx = barChartPlot.Parent;
          grid(plotAx, 'on');
          plotAx.Layer = 'top';
          
          title(titleText, 'Interpreter', 'none'); 
          xlabel('Jahr');
          ylabel(unit);
          legend(barChartPlot([2:5]), '0.05 - 0.25', '0.25 - 0.5', '0.5 - 0.75', '0.75 - 0.95', 'Location', 'Best');
          %set(fig,'PaperPositionMode','auto');
          set(fig, 'PaperUnits', 'inches');
          x_width=9.125;y_width=7.25;
          set(fig, 'PaperPosition', [0 0 x_width y_width]);
          saveas(fig, outJpgFileName, 'png');
        end
        
        function plotDensity(simulations, titleText, unit, scal, outJpgFileName)
          hor = size(simulations, 2);
          fig = figure;
          for t=1:hor
            subplot(ceil(hor/2), 2 , t);  
            histogram(scal * simulations(:, t), 'Normalization', 'probability');
            grid on
            title(['Jahr ', num2str(t)], 'Interpreter','none'); 
            xlabel(unit);
            ylabel('P');
          end
          %set(fig,'PaperPositionMode','auto');
          set(fig, 'PaperUnits', 'inches');
          suptitle(titleText);
          x_width=9.125;y_width=7.25;
          set(fig, 'PaperPosition', [0 0 x_width y_width]);
          saveas(fig, outJpgFileName, 'png');
        end
        
        function f = funRat(A, PN, V)
          f = (A - PN)./V;
        end  
        
        function s = sigma(q, mu, alpha)
          s = ((q - mu)/(sqrt(2) * erfinv(2 * alpha - 1)));
        end
        
        function plotPieChartStrategy(strat, labels, titleText, stratPlotJpgFileName)
            
          function pieChart0 = plotPieChartPortfolio(strat0)
              nonZeroPos = find(abs(strat0)>10^-3); 
              newStrat = strat0(nonZeroPos);
              newLabels = labels(nonZeroPos); 
              pieChart0 = pie(newStrat);
              h_legend = legend(newLabels, 'location', 'bestoutside');  
              set(h_legend, 'FontSize', 4);
          end
            
          hor = size(strat, 2);
          fig = figure;
          for t=1:hor
            subplot(ceil(hor/2), 2 , t);  
            plotPieChartPortfolio(strat(:,t));
            title(['Jahr ', num2str(t)], 'Interpreter','none'); 
          end
          set(fig, 'PaperUnits', 'inches');
          suptitle(titleText);
          x_width=9.125;y_width=7.25;
          set(fig, 'PaperPosition', [0 0 x_width y_width]);
          saveas(fig, stratPlotJpgFileName, 'png');
          
        end
        
        
                                                                      
    end    
    
    %public methods
    methods (Access = public)
    
        %constructor
        function obj = PensionFundProjection(pfa, pfl, co, pfr, pfrOpt, presDate)
            
            if not(listCompareQ(co.('Cat'), [getAssetsCategories(pfa); getPseudoLiabilitiesCategories(pfl)]))
               display('A&L:'); display([getAssetsCategories(pfa); getPseudoLiabilitiesCategories(pfl)]);
               display('CorrCat:'); display(co.('Cat'));
               msgID = 'PensionFundProjection:correlationIncompatibility';
               msg = strcat('Asset and liability categories not compatible with correlation category.');
               throw(MException(msgID,msg)); 
            end 
            
            if not(listCompareQ(getAssetsCategories(pfa), getAssetsCategories(pfr)))
               display('A:'); display(getAssetsCategories(pfa));
               display('Restr:'); display(getAssetsCategories(pfr));
               msgID = 'PensionFundProjection:restrictionsIncompatibility';
               msg = strcat('Asset categories not compatible with restriction category.');
               throw(MException(msgID,msg)); 
            end 
            
            if not(listCompareQ(getAssetsCategories(pfa), getAssetsCategories(pfrOpt)))
               display('A:'); display(getAssetsCategories(pfa));
               display('Restr Optimization:'); display(getAssetsCategories(pfrOpt));
               msgID = 'PensionFundProjection:restrictionsIncompatibility';
               msg = strcat('Asset categories not compatible with optimization restriction category.');
               throw(MException(msgID,msg)); 
            end 
           
            obj.pfass = pfa; 
            obj.pfliab = pfl;
            obj.corr = co;
            obj.pfrestr = pfr;
            obj.pfrestrOpt = pfrOpt;
            obj.presentDate = presDate;
            obj.simulTotRet = double.empty(0,0,0);
            obj.simulFunRat = double.empty(0,0);
            obj.assetSide = double.empty(0,0,0);
            obj.insuranceLiabilitySide = double.empty(0,0,0);
            obj.initialAssetStrategy = cellfun(@(arg) getAssetInitialMarketWeight(pfa, arg),getAssetsCategories(pfa));
            obj.assetStrategy = double.empty(0,0);
            obj.administrationCosts = double.empty(0);
            obj.inputsOpt = struct();
            obj.numSim = 5000;
            obj.numOpt = 30;
        end
        
        %Get/Set methods
        function  pfa = getPensionFundAssets(obj)
         pfa = obj.pfass;            
        end
        
        function  pfl = getPensionFundLiabilities(obj)
         pfl = obj.pfliab;            
        end
        
        function  pfr = getPensionFundRestrictions(obj)
         pfr = obj.pfrestr;            
        end
        
        function  pfrOpt = getPensionFundOptimizationRestrictions(obj)
         pfrOpt = obj.pfrestrOpt;            
        end
        
        function  presDate =  getPensionFundPresentDate(obj)
         presDate = obj.presentDate;            
        end
         
        function  numSim0 = getNumSim(obj)
         numSim0 = obj.numSim;            
        end
        
        function  setNumSim(obj, numSim0)
         obj.simulTotRet = double.empty(0,0,0);
	     clearCache(obj); 
	     obj.numSim = numSim0;            
        end
        
        function  numOpt0 = getNumberOptimalPortfolios(obj)
         numOpt0 = obj.numOpt;            
        end
        
        function  obj = setNumberOptimalPortfolios(obj, numOpt0) 
	     obj.numOpt = numOpt0;            
        end
        
        function  result = getCategories(obj)
         co = obj.corr;
         result = co.('Cat');
        end
        
        function  result = isCategoryPensionFundProjectionQ(obj, cat0)
         result = ismember(cat0, getCategories(obj));
        end
        
        function  result = getCorrelations(obj, t, cat1, cat2)
         co = obj.corr;   
         hor0 = length(setdiff(fieldnames(co),{'Cat'}));   
         if not(ismember(t,[1:hor0]))
               msgID = 'getCorrelations:badargException';
               msg = ['Time ', num2str(t), ' is not in interval 1,...,', num2str(hor0), '.'];
               throw(MException(msgID,msg)); 
         else 
            result = co.(['Period_', num2str(t)]);
            if nargin == 4
              if or(not(ismember(cat1,  getCategories(obj))), not(ismember(cat2,  getCategories(obj))))
                 msgID = 'getCorrelations:badargException';
                 msg = ['Either ', cat1, ' or ', cat2, ' is not in the asset/liability category list.'];
                throw(MException(msgID,msg)); 
              end
              result = result(catNum(obj, cat1), catNum(obj, cat2));
            end    
         end
        end
        
        function  result =  getTotalAssetInitialValue(obj)
         pfa = getPensionFundAssets(obj);
         result = sum(cellfun(@(arg) getAssetInitialMarketValue(pfa, arg), getAssetsCategories(pfa)));            
        end
        
        function  result =  getLHardCashFlowCosts(obj)
         pfl = getPensionFundLiabilities(obj);
         result = getLHardCfCosts(pfl);            
        end
        
        function  result =  getTotalNonInsuranceLiabilityInitialValue(obj)
         pfl = getPensionFundLiabilities(obj);
         result = sum(cellfun(@ (arg) getNonInsuranceLiabilityInitialValue(pfl, arg), getNonInsuranceLiabilitiesCategories(pfl)));          
        end
        
        function  result =  getTotalInsuranceRunOffCashFlows(obj)
         pfl = getPensionFundLiabilities(obj);
         result = getTotalInsuranceRunOffCashFlows(pfl);
        end
        
        function  result =  getTotalInsuranceCashFlows(obj)
         pfl = getPensionFundLiabilities(obj);
         result = getTotalInsuranceCashFlows(pfl);
        end
        
        function  result =  getVCashFlows(obj)
         %simulated insurance cash flows for every year (dummy
         %implementation while waiting for stochastic pension plan cash flows model)
         nS = getNumSim(obj);
         meanNetFlows = getTotalInsuranceCashFlows(obj);
         result = repmat(meanNetFlows, 2 * nS, 1);
        end
        
        function  result =  getExpectedTotalInsuranceCashFlows(obj, outputCsvFileName)
         result = mean(getVCashFlows(obj));      
         if nargin >=2
            hor = size(result,2);
            columnNamesValues = cell(1, hor);
            for t=1:hor
                columnNamesValues(t) = {['Year_', num2str(t)]};
            end  
              tableOutPut = obj.createGenericTableOutput(result, {'Insurance cash flow'}, columnNamesValues);  
              writetable(tableOutPut, outputCsvFileName, 'WriteRowNames', true, 'Delimiter', ';'); 
         end    
        end
        
        function  result =  getTotalInsuranceValues(obj)
         pfl = getPensionFundLiabilities(obj);
         result = sum(cellfunVectorOutput(@ (arg) getInsuranceLiabilityValues(pfl, arg), getInsuranceLiabilitiesCategories(pfl)));   
        end
        
        function  result =  getExpectedTotalNonInsuranceCashFlows(obj, outputCsvFileName)
         result = mean(getPNCashFlows(obj));
         if nargin >=2
            hor = size(result,2);
            columnNamesValues = cell(1, hor);
            for t=1:hor
                columnNamesValues(t) = {['Year_', num2str(t)]};
            end  
              tableOutPut = obj.createGenericTableOutput(result, {'Non-insurance cash flow'}, columnNamesValues);  
              writetable(tableOutPut, outputCsvFileName, 'WriteRowNames', true, 'Delimiter', ';'); 
         end    
        end
        
        function  result =  getExpectedTotalLiabilityCashFlows(obj)
         result = mean(getTotalLiabilityCashFlows(obj));    
        end
        
        function weights = getAssetStrategy(obj)
         weights = obj.assetStrategy;
        end
        
        function simulations = getSimulatedTotalReturns(obj)
            
         pfa = getPensionFundAssets(obj);
         pfl = getPensionFundLiabilities(obj);
         if isempty(obj.simulTotRet) 
            display('Simulation total return assets and liabilities for all periods...');
            assCat = getAssetsCategories(pfa);
            liabCat = getPseudoLiabilitiesCategories(pfl);
            cats = [assCat; liabCat];
            corrCats = getCategories(obj);
            cats2corrcatsInd = cellfun(@ (cat) find(ismember(corrCats, cat)), cats);
            retAss = cellfunVectorOutput(@ (cat) getAssetExpectedTotalReturns(pfa, cat), assCat);
            hor = length(retAss(1,:));
            retLiab = cellfunVectorOutput(@ (cat) getPseudoLiabilityExpectedTotalReturns(pfl, cat), liabCat); 
            volAss = cellfunVectorOutput(@ (cat) getAssetTotalReturnVolatilities(pfa, cat), assCat); 
            volLiab = cellfunVectorOutput(@ (cat) getPseudoLiabilityTotalReturnVolatilities(pfl, cat), liabCat);
            mus = [retAss; retLiab];
            vols = [volAss; volLiab];
            nS = getNumSim(obj);
            simulations = zeros(2*nS, hor, length(cats));
            rng(12345); %set random number geneterator seed for reproducibility
            for t=1:hor
                runningMu = mus(:,t);
                runningVols = vols(:,t);
                runningCorrs = getCorrelations(obj, t);
                runningCorrs = runningCorrs(cats2corrcatsInd, :);
                runningCorrs = runningCorrs(:, cats2corrcatsInd);
                runningSigma = diag(runningVols) * runningCorrs * diag(runningVols);
                [U, Lambda]= eig(runningSigma);
                %For highly correlated assets MATLAB can obtain light negative eigenvalues and complex eigenvectors...
                lambdas = diag(Lambda);
                lambdas(lambdas<10^-10.)=10^-10.;
                sqrtSigma = real(U*diag(sqrt(lambdas))*inv(U));
                runnigSim = mvnrnd(zeros(size(runningMu)), eye(size(runningSigma)), nS);
                runnigSim = [runnigSim;-runnigSim];  %variance reduction
                simulations(:, t, :) = runnigSim * sqrtSigma + repmat(runningMu',2*nS,1);
            end  
            obj.simulTotRet = simulations;
         else
            simulations = obj.simulTotRet;  
         end
         
        end
        
        function result = getTotalLiabilityCashFlows(obj)
          insuranceCF = getVCashFlows(obj); %dummy stochastic
          nonInsuranceCF = getPNCashFlows(obj); %stochastic
          pfl = getPensionFundLiabilities(obj); 
          totLCosts = getLAdministrationCosts(pfl) + getLHardCfCosts(pfl); %deterministic, total L costs including reorganization (Sanierung)
          result = repmat(totLCosts, size(insuranceCF,1),1) + insuranceCF + nonInsuranceCF;            
        end 
        
        function result = getTotalCashFlows(obj)
          insuranceCF = getVCashFlows(obj); %dummy stochastic
          nonInsuranceCF = getPNCashFlows(obj); %stochastic
          reOrg = getLHardCashFlowCosts(obj);%reorganization costs (Sanierung) 
          cost = getAdministrationCosts(obj);%administration costs A + L without reorganization (Sanierung) 
          result = repmat(reOrg+cost, size(insuranceCF,1),1) + insuranceCF + nonInsuranceCF;            
        end 
               
        function result = getAssetCashFlows(obj, outputCsvFileName)
           pfa = getPensionFundAssets(obj);
           assCat = getAssetsCategories(pfa);
           statValues = projectionAssetSide(obj);
           expectedAssetSide = [getTotalAssetInitialValue(obj), statValues(1, 1:end-1)];
           weightsA = getAssetStrategy(obj);
           runYields = cellfunVectorOutput(@ (arg) getAssetExpectedRunningYields(pfa, arg), assCat);
           result = repmat(expectedAssetSide, length(assCat), 1) .* weightsA .* runYields;
           
           if nargin == 2
              hor = size(statValues,2);
              columnNamesValues = cell(1, hor);
              for t=1:hor
                columnNamesValues(t) = {['Jahr_', num2str(t)]};
              end  
              text = ['ExpectedAssetCashFlows_', obj.assetStrategyName];
              tableOutPut = obj.createGenericTableOutput(result, assCat, columnNamesValues, text);  
              writetable(tableOutPut, outputCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
           end 
        end    
        
        function result = getTotalAssetCashFlows(obj)
           result = sum(getAssetCashFlows(obj));
        end    
        
        function statValues = projectionAssetSide(obj, cCsvFileName, aPlotJpgFileName, distPlotJpgFileName)

          cfs = getTotalCashFlows(obj);%includes all A + L costs and reorganization (Sanierung) costs
          hor = size(cfs, 2);
          incs = getAIncreaseFactor(obj);
          A0 = getTotalAssetInitialValue(obj);
          nS = getNumSim(obj);

          function as = proj(inc, cf)
            as = zeros(1, hor+1);
            as(1) = A0;
            for t=1:hor
             as(t+1) = as(t) * inc(t) - cf(t);
            end
            as = as(2:end);
          end  

          if isempty(obj.assetSide)
             aM = zeros(2*nS, hor); 
             for omega=1:2*nS
              aM(omega, :) = proj(incs(omega, :), cfs(omega, :));
             end
             obj.assetSide = aM;
          end   

          statValues = obj.getStatistics(obj.assetSide);
          text = ['MaketValuation_AssetSide_', obj.assetStrategyName];
          tableOutput = obj.createTableOutput(statValues, text);
          if nargin >= 2
             writetable(tableOutput, cCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Quantiles Asset Side: ', obj.assetStrategyName];
             unit = 'CHF';
             scal = 1.;
             obj.barPlot(statValues, titleText, unit, scal, aPlotJpgFileName);
          end  

          if nargin >=4
             titleText = ['Asset Side: ', obj.assetStrategyName]; 
             unit = 'CHF';
             scal = 1.;
             obj.plotDensity(obj.assetSide, titleText, unit, scal, distPlotJpgFileName);
          end


        end
        
        
        function statValues = projectionTotalReturnsAssetSide(obj, rCsvFileName, rPlotJpgFileName, distPlotJpgFileName)
          %asset only, no liabilities
           
          r = getAIncreaseFactor(obj) - 1;
          statValues = obj.getStatistics(r);
          text = ['Total_Returns_Assets_', obj.assetStrategyName];
          tableOutput = obj.createTableOutput(statValues, text);
          if nargin >= 2
             writetable(tableOutput, rCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Quantiles Total Returns Assets: ', obj.assetStrategyName];
             unit = '%';
             scal = 100.;
             obj.barPlot(statValues, titleText, unit, scal, rPlotJpgFileName);
          end  

          if nargin >=4
             titleText = ['Total Returns Assets: ', obj.assetStrategyName]; 
             unit = '%';
             scal = 100.; 
             obj.plotDensity(r, titleText, unit, scal, distPlotJpgFileName);
          end

        end
        
        function statValues = projectionCumulatedTotalReturnsAssetSide(obj, rCsvFileName, rPlotJpgFileName, distPlotJpgFileName)
          %asset only, no liabilities
           
          cumR = cumprod(getAIncreaseFactor(obj), 2) - 1;
          statValues = obj.getStatistics(cumR);
          text = ['Cumulated_Total_Returns_Assets_', obj.assetStrategyName];
          tableOutput = obj.createTableOutput(statValues, text);
          if nargin >= 2
             writetable(tableOutput, rCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Quantiles Cumulated Total Returns Assets: ', obj.assetStrategyName];
             unit = '%';
             scal = 100.;
             obj.barPlot(statValues, titleText, unit, scal, rPlotJpgFileName);
          end  

          if nargin >=4
             titleText = ['Cumulated Total Returns Assets: ', obj.assetStrategyName]; 
             unit = '%';
             scal = 100.;
             obj.plotDensity(cumR, titleText, unit, scal, distPlotJpgFileName);
          end

        end
        
        
        function statValues = projectionFundingRatio(obj, frCsvFileName, frPlotJpgFileName, distPlotJpgFileName)
                
          if isempty(obj.simulFunRat)
             projectionAssetSide(obj);
             pfl = getPensionFundLiabilities(obj);
             nILiabCat = getNonInsuranceLiabilitiesCategories(pfl);
             As = obj.assetSide;
             DKs = getTotalInsuranceValues(obj);
             DKs = DKs(2:end);
             PNs = sum(cellfunVectorOutput(@ (arg) getNonInsuranceLiabilityValuePrescribedEvolution(pfl, arg), nILiabCat));
             nS = getNumSim(obj);
             obj.simulFunRat = obj.funRat(As, repmat(PNs, 2*nS,1), repmat(DKs, 2*nS, 1)); 
          end   
            
          F = obj.simulFunRat;
          statValues = obj.getStatistics(F);
          text = ['FundingRatio_', obj.assetStrategyName];
          tableOutput = obj.createTableOutput(statValues, text);
          
          if nargin >= 2
             writetable(tableOutput, frCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Quantiles Funding Ratio: ', obj.assetStrategyName];
             unit = '%';
             scal = 100.;
             obj.barPlot(statValues, titleText, unit, scal, frPlotJpgFileName);
          end  

          if nargin >=4
             titleText = ['Funding Ratio: ', obj.assetStrategyName]; 
             unit = '%';
             scal = 100.;
             obj.plotDensity(F, titleText, unit, scal, distPlotJpgFileName);
          end

        end
        
        
        function statValues = projectionReportedInsuranceLiabilitySide(obj, vCsvFileName, vPlotJpgFileName, distPlotJpgFileName)
          %This method currently displays what has been read in from the XLSX Front End. There is no simulation for the book value of the
          %insurance liabilities and a normal distribution is assumed. To be improved with a stochastic model for the insurance liabilities
          
          pfl = getPensionFundLiabilities(obj);
          DKs = getTotalInsuranceValues(obj);
          DKs = DKs(2:end);
          quantDKs = getQuantilesReserves(pfl);
          quantDKs = flipud(quantDKs);
          alphas = quantDKs(:,1);
          sigmas = obj.sigma(quantDKs(1,2:end), DKs, alphas(1));
          statValues = [DKs; sigmas; quantDKs(:, 2:end)];
          text = 'Insurance_Liabilities';
          tableOutput = obj.createTableOutput(statValues, text);
          
          if nargin >= 2
             writetable(tableOutput, vCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = 'Quantiles Insurance-Liabilities';
             unit = 'CHF';
             scal = 1.;
             obj.barPlot(statValues, titleText, unit, scal, vPlotJpgFileName);
          end  

        end
        
        function statValues = projectionCashFlowInsurance(obj, cfCsvFileName, cfPlotJpgFileName, distPlotJpgFileName)
          %This method currently displays what has been read in from the XLSX Front End. There is no simulation for the cash flows of the
          %insurance liabilities and a normal distribution is assumed. To be improved with a stochastic model for the insurance liabilities
          pfl = getPensionFundLiabilities(obj);
          cf = getExpectedTotalInsuranceCashFlows(obj);
          quantCfs = getQuantilesInsuranceCashFlows(pfl);
          %Inversion of cash flow signs
          cf = -cf;
          quantCfs(:,1) = 1 - quantCfs(:,1);
          quantCfs(:,2:end) = -quantCfs(:,2:end);
          alphas = quantCfs(:,1);            
          sigmas = obj.sigma(quantCfs(1,2:end), cf, alphas(1));
          statValues = [cf; sigmas; quantCfs(:, 2:end)];
          text = 'Cash_Flows_Insurance_Liabilities';
          tableOutput = obj.createTableOutput(statValues, text);
          
          if nargin >= 2
             writetable(tableOutput, cfCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = 'Quantiles Cash Flows Insurance-Liabilities';
             unit = 'CHF';
             scal = 1.;
             obj.barPlot(statValues, titleText, unit, scal, cfPlotJpgFileName);
          end  

        end
                
        
        function f = getInitialFundingRatio(obj, outputCsvFileName) 
           DKs = getTotalInsuranceValues(obj); 
           f = obj.funRat(getTotalAssetInitialValue(obj), getTotalNonInsuranceLiabilityInitialValue(obj), DKs(1));
           if nargin >= 2 
              columnNamesValues = {'Year_0'};
              rowLabels = {'Funding Ratio'};
              tableOutPut = obj.createGenericTableOutput(f, rowLabels, columnNamesValues);  
              writetable(tableOutPut, outputCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');    
           end
        end   
        
        function result = getTargetValueFluctuationReserves(obj, alpha, targetVFRCsvFileName)
          
            if isempty(obj.simulFunRat)  
             projectionFundingRatio(obj);
            end 
           
           DKs = getTotalInsuranceValues(obj);
           cfs = getExpectedTotalLiabilityCashFlows(obj);
           RVs = (DKs(2:end)+ cfs)./DKs(1:end-1) - 1;
           retAs = getAIncreaseFactor(obj) - 1;
           VaR = quantile(retAs, alpha, 1);
           CVaR = [];
           for t = 1:size(retAs,2)
            CVaR = [CVaR, 1/alpha * mean(retAs(:,t) .* (retAs(:,t)<=quantile(retAs(:,t), alpha)))];
           end
           result = [- VaR + RVs; (1 + RVs) ./ (1 + VaR) - 1; (1 + RVs) ./ (1 + CVaR) - 1];
           
           if nargin >= 3
              hor = size(result,2);
              columnNamesValues = cell(1, hor);
              for t=1:hor
                columnNamesValues(t) = {['Year_', num2str(t)]};
              end  
              rowLabels = {['Proxy-VaR-Target-VFR ', num2str(100*(1-alpha)), '%'], ['VaR-Target-VFR ', num2str(100*(1-alpha)), '%'], ['CVaR-Target-VFR ', num2str(100*(1-alpha)), '%']};
              text = obj.assetStrategyName;
              tableOutPut = obj.createGenericTableOutput(result, rowLabels, columnNamesValues, text);  
              writetable(tableOutPut, targetVFRCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');    
           end
        
        end 
        
         function result = getCumulatedTargetValueFluctuationReserves(obj, alpha, targetVFRCsvFileName)
          
            if isempty(obj.simulFunRat)  
             projectionFundingRatio(obj);
            end 
           
           DKs = getTotalInsuranceValues(obj);
           cfs = getExpectedTotalLiabilityCashFlows(obj);
           RVs = cumprod((DKs(2:end)+ cfs)./DKs(1:end-1), 2) - 1;
           retAs = cumprod(getAIncreaseFactor(obj), 2) - 1;
           VaR = quantile(retAs, alpha, 1);
           CVaR = [];
           for t = 1:size(retAs,2)
            CVaR = [CVaR, 1/alpha * mean(retAs(:,t) .* (retAs(:,t)<=quantile(retAs(:,t), alpha)))];
           end
           result = [- VaR + RVs; (1 + RVs) ./ (1 + VaR) - 1; (1 + RVs) ./ (1 + CVaR) - 1];
           
           if nargin >= 3
              hor = size(result,2);
              columnNamesValues = cell(1, hor);
              for t=1:hor
                columnNamesValues(t) = {['Year_', num2str(t)]};
              end  
              rowLabels = {['Proxy-VaR-Target-VFR ', num2str(100*(1-alpha)), '%'], ['VaR-Target-VFR ', num2str(100*(1-alpha)), '%'], ['CVaR-Target-VFR ', num2str(100*(1-alpha)), '%']};
              text = obj.assetStrategyName;
              tableOutPut = obj.createGenericTableOutput(result, rowLabels, columnNamesValues, text);  
              writetable(tableOutPut, targetVFRCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');    
           end
        
        end 
         
        function terminalWeightsA = checkRestrictionsOverTime(obj)
          pfa = getPensionFundAssets(obj);
          assCat = getAssetsCategories(pfa);
          rets =  cellfunVectorOutput(@ (arg) getAssetExpectedTotalReturns(pfa, arg), assCat);
          weightsA = getAssetReportedStrategy(obj);
          hor = size(weightsA, 2);
          terminalWeightsA = zeros(size(weightsA));
          for t=1:hor
           terminalWeightsA(:,t) = relW(weightsA(:,t) .* (1 + rets(:,t)));
          end
          restrFn = getRestrictionsFunction(obj);
          test = @ (year) all(restrFn(terminalWeightsA(:, year))<=0);
          for t=1:hor
              if not(test(t))
                 display(['RESTRICTIONS NOT FULLFILLED BETWEEN YEAR ', num2str(t - 1), ' AND ', num2str(t), '!']);
              else
                 display(['Restrictions fullfilled between year ', num2str(t - 1), ' und ', num2str(t), '.']);
              end    
          end
        end 
        
        function alphas = getReportedFundingRatioProbability(obj, F, alphasCsvFileName)
            
          if isempty(obj.simulFunRat)  
             projectionFundingRatio(obj)
          end
          
          f = obj.simulFunRat;
          alphas = mean(f<=F);
          
          if nargin == 3
             hor = size(alphas,2);
             columnNamesValues = cell(1, hor);
             for t=1:hor
              columnNamesValues(t) = {['Year_', num2str(t)]};
             end  
             rowNamesValues = {['P[f_t<=', num2str(100*F),'%]']};
             text = obj.assetStrategyName;
             tableOutput = obj.createGenericTableOutput(alphas, rowNamesValues, columnNamesValues, text); 
             writetable(tableOutput, alphasCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 
          
        end 
        
        function alphas = projectionDefaultProbabilities(obj, fLevels, alphasCsvFileName, alphasPlotJpgFileName)
         
          alphas = [];
          for level=1:length(fLevels)
           alphas = [alphas; getReportedFundingRatioProbability(obj, fLevels(level))];
          end
 
          if nargin >= 3
             hor = size(alphas,2);
             columnNamesValues = cell(1, hor);
             for t=1:hor
              columnNamesValues(t) = {['Year_', num2str(t)]};
             end  
             rowNamesValues = {};
             for level=1:length(fLevels)
              rowNamesValues = [rowNamesValues; {['P[f_t<=', num2str(100 * fLevels(level)),'%]']}];
             end
             text = obj.assetStrategyName;
             tableOutput = obj.createGenericTableOutput(alphas, rowNamesValues, columnNamesValues, text); 
             writetable(tableOutput, alphasCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=4
             titleText = ['Probability that the funding ratio falls below a threshold: ', text];
             fig = figure;
             plotAlphas = plot(100 * alphas', 'LineWidth', 1.25);
             grid on
             title(titleText, 'Interpreter','none'); 
             xlabel('Jahr');
             ylabel('P[%]');
             labels = {};
             for level=1:length(fLevels)
              labels = [labels; {[num2str(100*fLevels(level)),'%']}];
             end
             legend(labels, 'Location', 'Best');
             %set(fig,'PaperPositionMode','auto');
             set(fig, 'PaperUnits', 'inches');
             x_width=9.125;y_width=7.25;
             set(fig, 'PaperPosition', [0 0 x_width y_width]);
             saveas(fig, alphasPlotJpgFileName, 'png');
          end  

        end
        
        
        function bs = getBalanceSheetProjection(obj, bsCsvFileName, bsPlotJpgFileName)
          %Regulatory balance sheet projection at expected value level  
          
          pfl = getPensionFundLiabilities(obj);
          nILiabCat = getNonInsuranceLiabilitiesCategories(pfl);
          a0 = getTotalAssetInitialValue(obj);
          pn0 = getTotalNonInsuranceLiabilityInitialValue(obj);
          DK = getTotalInsuranceValues(obj);
          pn = sum(cellfunVectorOutput(@ (arg) getNonInsuranceLiabilityValuePrescribedEvolution(pfl,arg), nILiabCat));
          pn = [pn0, pn];          
          if isempty(obj.assetSide)
             projectionReportedAssetSide(obj);
          end
          a = mean(obj.assetSide);
          a = [a0, a];
          s = a - pn - DK;
          bs = [a;DK;pn;s];
          
          defaultTimes = find(s<=0);
          if not(isempty(defaultTimes))
              display(['Warning! Default in year ', num2str(defaultTimes(1)), ' expected!']);
          end
          if nargin >= 2
             hor = size(bs,2);
             columnNamesValues = cell(1, hor);
             for t=1:hor
              columnNamesValues(t) = {['Jahr_', num2str(t-1)]};
             end  
             rowNamesValues = {'Assets'; 'Insurance Liabilities'; 'Non-Insurance Liabilities'; 'Value Fluctuation Reserve'};
             text = obj.assetStrategyName;
             tableOutput = obj.createGenericTableOutput(bs, rowNamesValues, columnNamesValues, text); 
             writetable(tableOutput, bsCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Expected Balance Sheet: ', text];
             unit = 'CHF';
             %fig = figure;
             groupLabels = {}; 
             for t=1:hor
                 groupLabels = [groupLabels, {num2str(t-1)}];
             end
             stackData = zeros(hor, 2, 3);
             for t=1:hor
                stackData(t,1,1)= bs(1, t);%a
                %rearrange order of DK, pn, s
                for kk=2:4
                 stackData(t,2,kk)= bs(3-(kk-3),t);
                end 
             end    
             %rearrange order labels
             labels = [rowNamesValues(1);flipud(rowNamesValues(2:end))];
             barChartPlot= plotBarStackGroups(stackData, groupLabels);
             %external function plotBarStackGroups created its own figure
             %already
             fig = gcf;
             % have to change colors of chart and legend at the same time
             %display(barChartPlot.FaceColor)
             barChartPlot(1,1).FaceColor = [255,0,0]/255;
             barChartPlot(2,1).FaceColor = [255,0,0]/255;
             barChartPlot(1,2).FaceColor = [0,176,240]/255;
             barChartPlot(2,2).FaceColor = [0,176,240]/255;
             barChartPlot(1,3).FaceColor = [255,255,0]/255;
             barChartPlot(2,3).FaceColor = [255,255,0]/255;
             barChartPlot(1,4).FaceColor = [0,176,80]/255; 
             barChartPlot(2,4).FaceColor = [0,176,80]/255; 
%              bl = barChartPlot.BaseLine;
%              bl.Visible = 'off';
             grid on
             title(titleText, 'Interpreter','none'); 
             xlabel('Year');
             ylabel(unit);
             legend(labels, 'Location', 'Best');
             set(fig, 'PaperUnits', 'inches');
             x_width=9.125;y_width=7.25;
             set(fig, 'PaperPosition', [0 0 x_width y_width]);
             saveas(fig, bsPlotJpgFileName, 'png');
          end  

        end
        
        function strats = getPieChartAssetStrategy(obj, stratCsvFileName, stratPlotJpgFileName)
          
          pfa = getPensionFundAssets(obj);
          strats = getAssetStrategy(obj);
          labels = getAssetsCategories(pfa);
          
          if nargin >= 2
             hor = size(strats,2);
             columnNamesValues = cell(1, hor);
             for t=1:hor
              columnNamesValues(t) = {['Year_', num2str(t)]};
             end  
             rowNamesValues = labels;
             text = obj.assetStrategyName;
             tableOutput = obj.createGenericTableOutput(strats, rowNamesValues, columnNamesValues, text); 
             writetable(tableOutput, stratCsvFileName, 'WriteRowNames', true, 'Delimiter', ';');
          end 

          if nargin >=3
             titleText = ['Asset Allocation ', text];
             obj.plotPieChartStrategy(strats, labels, titleText, stratPlotJpgFileName);
          end  

        end
        
        
        function calcAssetMarketStrategy(obj, altStrat, altStratName)
          
          if nargin <= 2
             altStratName = 'Alternative';
          end
          obj.assetStrategyName = altStratName;
          pfa = getPensionFundAssets(obj);
          clearCache(obj);
          assCat = getAssetsCategories(pfa);
          nAss = length(assCat);
          candWeights = altStrat(1:end-1,1);
          restrFn =  getRestrictionsFunction(obj);
          
          test = (abs(sum(candWeights)-1)< 10^-4)&& all(restrFn(candWeights)<=0);
          
          if not(test)
             msgID = 'calcAssetMarketStrategy:badargException';
             msg = ['Alternative allocation at time t=0 does not satisfy restrictions.'];
             display(msg);
             %throw(MException(msgID,msg)); 
          end
          
          obj.initialAssetStrategy = candWeights;
          obj.initialAssetStrategyR = candWeights;
          %New initial market and reported asset strategies are equal because of rebalancing
          candFMV = altStrat(1:end-1, 2:end);
          for j=1:nAss
           setAssetPrescribedFutureMarketValues(pfa, assCat{j}, candFMV(j,:));
          end
          candRETC = altStrat(end, 2:end);
          setAHardCfCosts(pfa, candRETC);
          initialize(obj, 1);           
          
        end    
        
       
% TO BE IMPLEMENTED 
%         getOptimizationInequalityRestrictions;
%         getOptimalTacticalBounds;
%         getALMOptimalAllocation;
%         getALMEfficientFrontier;

% NOT IMPLEMENTED BECAUSE TILL NOW NOT NEEDED
%         getMarketBasedFundingRatioProbability;probably not needed
%         getPieChartReportedInsuranceLiability;
%         getPieChartReportedNonInsuranceLiability;

        
    end
    
    %private methods
    methods (Access = public) %Remember to set this to private when the implementation is over
        
        function clearCache(obj)
            
            obj.simulTotRet = double.empty(0,0,0);
            obj.simulFunRat = double.empty(0,0);
            obj.assetSide = double.empty(0,0,0);
            obj.insuranceLiabilitySide = double.empty(0,0,0);
            obj.assetStrategy = double.empty(0,0);
            obj.administrationCosts = double.empty(0);
            pfa = getPensionFundAssets(obj);
            obj.initialAssetStrategy = cellfun(@(arg) getAssetInitialMarketWeight(pfa, arg),getAssetsCategories(pfa));
              
        end    
        
        function catNum0 = catNum(obj, cat0)
            catNum0 = find(ismember(getCategories(obj), cat0));    
        end
        
        function initialize(obj, testFeasibility)
            
            pfl = getPensionFundLiabilities(obj);
            pfa = getPensionFundAssets(obj); 
            A0 = getTotalAssetInitialValue(obj);
             
            display('Initialization...');  
            cf = mean(getVCashFlows(obj) + getPNCashFlows(obj));         
            hor = length(cf);
            assCat = getAssetsCategories(pfa);
            for t=1:hor
             assPresFMV.(['Time_',num2str(t)]) = getAssetsWithPrescribedFutureValuesCategories(pfa, t);
             assNotPresFMV.(['Time_',num2str(t)]) = getAssetsWithoutPrescribedFutureValuesCategories(pfa, t);
            end
            %why using expected and not simulated asset returns to compute costs?
            rets = cellfunVectorOutput(@ (cat) getAssetExpectedTotalReturns(pfa,cat), assCat);
            
            function out = fmvs(time)
              function res = getFMV(arg)
                 vec = setAssetPrescribedFutureValues(pfa, arg);
                 res = vec(time);
              end
              out = cellfun(@getFMV, assPresFMV.(['Time_',num2str(time)])); 
            end
            
            rOrg =  getLHardCashFlowCosts(obj); %reorganization cashflows (Sanierung)
            reCosts = getAHardCfCosts(pfa);% assets strategy costs, hard cash flows,  typically for real estate
            
            lCosts = getLAdministrationCosts(pfl); %pension fund administration costs WITHOUT reorganization cashflows (Sanierung)
            aCostsPer = getAAdministrationCosts(pfa);
            weightsM0 = obj.initialAssetStrategy;
            weightsM = zeros(length(weightsM0), hor); 
            weightsM(:,1) = weightsM0;
             
            posPresFMV = @(t) find(ismember(assCat, assPresFMV.(['Time_',num2str(t)])));
            posNotPresFMV = @(t) find(ismember(assCat, assNotPresFMV.(['Time_',num2str(t)])));
             
            %relWeigthsR0 = relW(weightsR0(posNotPresFMV(1)));
            %relWeigthsM0 = relW(weightsM0(posNotPresFMV(1)));
            posSec = find(ismember(assCat,  getSecuritiesCategories(pfa)));
            secWeight0 = sum(weightsM0(posSec));
            secWeight = zeros(1, hor);
            secWeight(1) = secWeight0;

            incA0 = 1 + weightsM0'*rets(:,1);
            incA = zeros(1, hor);
            incA(1) = incA0;
            As = zeros(1, hor+1);
            As(1) = A0;
            aCosts = zeros(1, hor);
            aCosts(1) = aCostsPer(1) * secWeight0 * A0 + reCosts(1);
            for t=1:(hor-1)
              As(t + 1) = As(t) * incA(t) - cf(t) - lCosts(t) - aCosts(t) - rOrg(t); 
              secWeight(t + 1) = sum(weightsM(posSec, t + 1));
              aCosts(t + 1) = aCostsPer(t + 1) * secWeight(t + 1) * As(t + 1) + reCosts(t + 1);
              weightsM(posNotPresFMV(t), t + 1) = relW(weightsM(posNotPresFMV(t), t)) * (As(t + 1) - sum(fmvs(t)))/As(t + 1);
              weightsM(posPresFMV(t), t+1)= fmvs(t)/As(t + 1);
              incA(t + 1) = 1 + weightsM(:, t + 1)'*rets(:, t + 1);
            end
            
            obj.administrationCosts = aCosts + lCosts;
            obj.assetStrategy = weightsM;
            
            
            if nargin == 2
               restrFn =  getRestrictionsFunction(obj);
               test = @ (year) all(restrFn(weightsM(:, year))<=0);
               if testFeasibility
                  for t=1:hor 
                    if not(test(t))
                       display(['Unfeasible assetStrategy at time ', num2str(t - 1), ':']);
                       weightsM(:, t)
                       restrFn(weightsM(:, t))<=0
                    end
                  end
               end
            end
               
        end    
        
        function result = getRestrictionsFunction(obj)
            
            pfr =  getPensionFundRestrictions(obj);
            pfa = getPensionFundAssets(obj);
            assCat = getAssetsCategories(pfa);
            resCat = getAssetsCategories(pfr);
            mixedBoundsMatrix = getNonElementaryConstraints(pfr);
            mixedUpperBounds = mixedBoundsMatrix(end,:);
            mixedCoeffs = mixedBoundsMatrix(1:end-1,:);
            permutedMixedCoeff = mixedCoeffs(cellfun(@ (cat) find(ismember(cat, resCat)), assCat), :);
            
            function out = fn(syms)
             lowerBounds = cellfun(@ (cat) getLowerBound(pfr, cat), assCat) - syms;
             upperBounds = syms -  cellfun(@ (cat) getUpperBound(pfr, cat), assCat);
             mixedBounds = permutedMixedCoeff'*syms - mixedUpperBounds';
             out = [lowerBounds;upperBounds; mixedBounds];
            end
            
            result = @fn;
        
        end
        
        function result = getAReturns(obj)
          pfa = getPensionFundAssets(obj);  
          assCat = getAssetsCategories(pfa);
          assInd= cellfun(@ (arg) catNum(obj, arg), assCat);
          simulations = getSimulatedTotalReturns(obj);
          result = simulations(:, :, assInd);
        end   
        
        function result = getPNReturns(obj)
          pfl = getPensionFundLiabilities(obj);  
          nilCat = getNonInsuranceLiabilitiesCategories(pfl);
          nilInd= cellfun(@ (arg) catNum(obj, arg), nilCat);
          simulations = getSimulatedTotalReturns(obj);
          result = simulations(:, :, nilInd);
        end  
        
        function result = getAIncreaseFactor(obj)
          weightsM = getAssetStrategy(obj);
          simRetA = getAReturns(obj);
          T = size(simRetA,2);
          nS = getNumSim(obj); 
          for t=1:T
              for omega = 1:2*nS
               result(omega, t) = 1 + squeeze(simRetA(omega, t, :))' * weightsM(:,t);
              end
          end
        end  
        
        function result = getPNCashFlows(obj)
          pfl = getPensionFundLiabilities(obj);
          cfRelNILiabCat = getNonInsuranceCashFlowRelevantLiabilitiesCategories(pfl);
          redInc = getCfRelevantPNIncreaseFactor(obj);
          nilME = cellfunVectorOutput(@ (arg) getNonInsuranceLiabilityValuePrescribedEvolution(pfl, arg), cfRelNILiabCat);
          nS = getNumSim(obj); 
          
          if isempty(cfRelNILiabCat)
              redPN = zeros(1, size(nilME,2)+1);
          else    
              redPN0 = sum(cellfun(@ (arg) getNonInsuranceLiabilityInitialValue(pfl, arg), cfRelNILiabCat));
              redPN = sum(cellfunVectorOutput(@ (arg) getNonInsuranceLiabilityValuePrescribedEvolution(pfl, arg), cfRelNILiabCat));
              redPN = [redPN0, redPN];
          end
          result = redInc .* repmat(redPN(1:end-1), 2*nS, 1) - repmat(redPN(2:end),  2*nS, 1);  
        end
        
        function result = getCfRelevantPNIncreaseFactor(obj)
          pfl = getPensionFundLiabilities(obj);
          nS = getNumSim(obj); 
          cfRelNILiabCat = getNonInsuranceCashFlowRelevantLiabilitiesCategories(pfl);
          nILiabCat = getNonInsuranceLiabilitiesCategories(pfl);
          weightsPN0 = cellfun(@(arg) getNonInsuranceLiabilityInitialWeight(pfl, arg), nILiabCat);
          matrixWeigthsPN = cellfunVectorOutput(@ (arg) getNonInsuranceLiabilityWeightPrescribedEvolution(pfl, arg), nILiabCat);
          matrixWeigthsPN = [weightsPN0, matrixWeigthsPN(:, 1:end-1)];
          zeroOne = cellfun(@ (arg) ismember(arg, cfRelNILiabCat), nILiabCat);
          redMatrixWeigthsPN = repmat(zeroOne, 1, size(matrixWeigthsPN,2)) .* matrixWeigthsPN;
          
          T = size(redMatrixWeigthsPN, 2);
          for t=1:T
              redMatrixWeigthsPN(:, t) = relW(redMatrixWeigthsPN(:, t));
          end 
          pnRets = getPNReturns(obj);
          result = zeros(2*nS, T);
          for omega=1:2*nS
             for t=1:T
              result(omega, t) = 1 + redMatrixWeigthsPN(:,t)' * squeeze(pnRets(omega, t, :));
             end 
          end
        end    
        
        function result = getAdministrationCosts(obj)  
          if isempty(obj.administrationCosts)
             initialize(obj, logical(1));
          end
          result = obj.administrationCosts;
        end 
   
        function result = getPNIncreaseFactor(obj)
          pfl = getPensionFundLiabilities(obj);
          nS = getNumSim(obj); 
          nILiabCat = getNonInsuranceLiabilitiesCategories(pfl);
          weightsPN0 = cellfun(@ (arg) getNonInsuranceLiabilityInitialWeight(pfl, arg), nILiabCat);
          matrixWeigthsPN = cellfunVectorOutput(@(arg) getNonInsuranceLiabilityWeightPrescribedEvolution(pfl, arg), nILiabCat);
          matrixWeigthsPN = [weightsPN0, matrixWeigthsPN(:,1:end-1)];
          T = size(matrixWeigthsPN, 2);
          pnRets = getPNReturns(obj);
          result = zeros(2*nS, T);
          for omega=1:2*nS
             for t=1:T
              result(omega, t) = 1 + matrixWeigthsPN(:,t)' * squeeze(pnRets(omega, t, :));
             end 
          end
        end
        
          
    end
    
end

