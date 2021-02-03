function strPerf = evalPredictions(predicted, GT, flagDisplay, perfs)
    strPerf = evalsinglePred(); strPerf.mean = NaN;
    ff = fieldnames(strPerf);
    
    if nargin > 1
        if nargin<4; perfs = {'auc'}; end
        if nargin<3; flagDisplay = false; end
        numclasses = numel(unique(GT));
        if (numclasses==2)
            if numel(predicted) == numel(GT)
               strPerf = evalsinglePred( predicted, GT);
            else; numRounds = size(predicted,2); 

                for nr = 1: numRounds
                    strR = evalsinglePred(predicted(:,nr), GT);
                    for f = 1:numel(ff); if ~isa(strPerf.(ff{f}), 'cell')
                       strPerf.(ff{f}) = ((nr-1)*strPerf.(ff{f})+strR.(ff{f}))/nr; end; end
                end
            end
        else
            % se ci sono più classi uso ogni classe come se fosse la
            % positiva e le altre come se fossero le negative. Poi faccio
            % la media delle performance
            for class = 1:numclasses
                  GTNow = GT ==  class;
                  predNow = predicted == class;
                  strPerfMean = strPerf;
                    if numel(predNow) == numel(GTNow)
                       strPerf = evalsinglePred(predNow, GTNow);
                    else; numRounds = size(predNow,2); 

                        for nr = 1: numRounds
                            strR = evalsinglePred(predNow(:,nr), GTNow);
                            for f = 1:numel(ff); if ~isa(strPerf.(ff{f}), 'cell')
                               strPerf.(ff{f}) = ((nr-1)*strPerf.(ff{f})+strR.(ff{f}))/nr; end; end
                        end
                    end
                  for f = 1:numel(ff); if ~isa(strPerf.(ff{f}), 'cell')
                      strPerfMean.(ff{f}) = strPerfMean.(ff{f})+strPerf.(ff{f}); end
                  end 
                    
            end
            for f = 1:numel(ff); if ~isa(strPerf.(ff{f}), 'cell')
                strPerfMean.(ff{f}) = strPerfMean.(ff{f})/numclasses; end
            end 
        end

        if  flagDisplay
            disp(['sens = ' num2str(strPerf.sens)  ...
                  ', spec  = ' num2str(strPerf.spec) ...
                  ', Acc = ' num2str(strPerf.acc)   ...
                  ', PPV = ' num2str(strPerf.ppv)   ...
                  ', NPV = ' num2str(strPerf.npv)   ...
                  ', F1  = ' num2str(strPerf.f1)   newline ...
                  'Critical Success Index/Threat score (TS) = ' num2str(strPerf.ts)   newline ...
                  'Informedness/Bookmaker Informedness (BM) = ' num2str(strPerf.bm)   ]);
        end  
        
        meanV = 0; 
        for np = 1: numel(perfs)
            if isfinite(strPerf.(perfs{np})); meanV = meanV + strPerf.(perfs{np}); 
              end; end
        strPerf.mean =  meanV;
    end
    
end
      
      
