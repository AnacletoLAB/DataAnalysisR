function strPerf = evalsinglePred(predicted, GT)
    
    strPerf.perfmeas ={'sens', 'spec', 'acc', 'PPV', 'NPV', 'F1', 'TS', 'BM', 'AUC', 'FP', 'FN'};
    strPerf.all = zeros(1,numel(strPerf.perfmeas));
    strPerf.sens = 0;
    strPerf.spec = 0;
    strPerf.acc = 0;
    strPerf.ppv = 0;
    strPerf.npv = 0;
    strPerf.f1 = 0;
    strPerf.ts = 0;
    strPerf.bm = 0;
    strPerf.auc = 0;
    strPerf.FP = 0;
    strPerf.FN = 0;
    
    if nargin == 2
        precdigit = 2;
        indDel = isnan(predicted) | isnan(GT);
        predicted(indDel) = []; GT(indDel) = [];
        crossCheck = confusionmat(categorical(round(predicted)), categorical(double(GT)));
        %[TP, FP, TN, FN] = eval(labelsBefore, labels)
        P = sum(crossCheck(:,2));
        N = sum(crossCheck(:,1));

        TP = crossCheck(2,2); sens = round(TP/P, precdigit); % sensitivity/recall/TruePositiveRate
        TN = crossCheck(1,1); spec = round(TN/N, precdigit); % specificity/selectivity/TrueNegativeRate
        FP = crossCheck(2,1); FPR = round(FP/N, precdigit); % fall-out/false positive rate
        FN = crossCheck(1,2); FNR = round(FN/P, precdigit); %miss rate/false negative rate

        PPV = round(TP/(TP+FP), precdigit); % precision/positive predictive value
        NPV = round(TN/(TN+FN), precdigit); % negative predictive value

        TS = round(TP/(TP+FN+FP), precdigit); % Critical Success Index/Threat score

        ACC = round((TP+TN)/(P+N), precdigit); % accuracy

        F1 = round((2*TP)/(2*TP+FP+FN), precdigit); % F1 score (harmonic mean of precision and sensitivity)

        BM = sens + spec -1; %Informedness or Bookmaker Informedness (BM)

        [~,~,~,auc]= perfcurve(GT,double(predicted),1);
        auc = round(auc, precdigit);
        strPerf.all = [sens, spec, ACC, PPV, NPV, F1, TS, BM, auc, FP, FN];
        strPerf.sens = sens;
        strPerf.spec = spec;
        strPerf.f1 = F1;
        strPerf.acc = ACC;
        strPerf.ppv = PPV;
        strPerf.npv = NPV;
        strPerf.ts = TS;
        strPerf.bm = BM;
        strPerf.auc = auc;
        strPerf.FP= FP;
        strPerf.FN = FN;
    end


end