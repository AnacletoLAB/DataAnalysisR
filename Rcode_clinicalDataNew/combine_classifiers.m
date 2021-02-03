function res = combine_classifiers()
    rng(1);
    data = readtable(['..' filesep 'data_covnet_score_selected-imputed_501.csv'], 'ReadVariableNames' , true);
    LABEL = data.LABEL_2;
    data.LABEL_2 = [];
    folds = table2array(readtable(['..' filesep 'folds.txt'], 'ReadVariableNames' , false))+1;
    data_norm = (table2array(data)-repmat(min(table2array(data)), size(data,1),1));
    data_norm = data_norm./repmat(max(data_norm), size(data_norm,1),1);
    idx_cat = 1:6;
    idx_sat = 7:9;
    idx_blood_test = 10:17;
    idx_radio = 18:21;
    preds_labels = [];
    str_optimizer = 'bayesopt';
    cost = [0 0.75; ... 
                0.35 0];
    for numf = 1:size(folds,1)
       idx_test = folds(numf,:);
       data_test = data_norm(idx_test,:);
       lab_test = LABEL(idx_test);
       
       data_train = data_norm;
       data_train(idx_test,:) = [];
       lab_train = LABEL;
       lab_train(idx_test,:) = [];
      
       idx_tree_train = 1:size(data_train,1);
       
       idx_svm_train = 1:size(data_train,1);
       idx_svm_train(idx_tree_train) = [];
       
       data_tree_train = data_train(idx_tree_train, :);
       lab_tree_train = lab_train(idx_tree_train, :);
       
       data_svm_train = data_train(idx_svm_train, :);
       lab_svm_train = lab_train(idx_svm_train, :);
       
       
       cv_bayes = cvpartition(lab_svm_train,'KFold',10);
       
       
        struct_opt = struct('Optimizer',str_optimizer,'MaxObjectiveEvaluations', 100, 'Useparallel', true, ...
           'cvpartition', cv_bayes);
       best_svm_cat = fitcsvm(data_svm_train(:, idx_cat), lab_svm_train, ...
           'Cost', cost, ...
           'OptimizeHyperparameters', 'all', ...
           'HyperparameterOptimizationOptions',struct_opt );

       best_svm_radio = fitcsvm(data_svm_train(:, idx_radio), lab_svm_train, ...
           'Cost', cost, ...
           'OptimizeHyperparameters', 'all', ...
           'HyperparameterOptimizationOptions', struct_opt);

       best_svm_blood = fitcsvm(data_svm_train(:, idx_blood_test), lab_svm_train,...
           'Cost', cost, ...
           'OptimizeHyperparameters', 'all', ...
           'HyperparameterOptimizationOptions', struct_opt);

       
       best_svm_sat = fitcsvm(data_svm_train(:, idx_sat), lab_svm_train, ...
           'Cost', cost, ...
           'OptimizeHyperparameters', 'all', ...
           'HyperparameterOptimizationOptions', struct_opt);

       [labs_cat_t, score_cat_t] = predict(best_svm_cat, data_tree_train(:, idx_cat));
       [labs_radio_t, score_radio_t] =  predict(best_svm_radio, data_tree_train(:, idx_radio));
       [labs_blood_t, score_blood_t] = predict(best_svm_blood,data_tree_train(:, idx_blood_test));
       [labs_sat_t, score_sat_t] = predict(best_svm_sat, data_tree_train(:, idx_sat));
       
       labels_pred_svm = [fix(score_cat_t) fix(score_radio_t) fix(score_blood_t) fix(score_sat_t)];
         
       tree_prediction = fitctree(labels_pred_svm, lab_tree_train);
       
       [labs_cat, score_cat] = predict(best_svm_cat, data_test(:, idx_cat));
       [labs_radio, score_radio] =  predict(best_svm_radio, data_test(:, idx_radio));
       [labs_blood, score_blood] = predict(best_svm_blood, data_test(:, idx_blood_test));
       [labs_sat, score_sat] = predict(best_svm_sat, data_test(:, idx_sat));
       
       
       labels_pred_test = [fix(score_cat) fix(score_radio) fix(score_blood) fix(score_sat)];
         
       preds  = predict(tree_prediction, labels_pred_test);
       preds_labels = [preds_labels;  preds lab_test];
       
       close all
       
    end
    %[SVMModel,ScoreParameters] = fitPosterior(SVMModel);
    
    evalPredictions(preds_labels(:,1), preds_labels(:,2))
end



function score = fix(score)
    score(score<-1) = -1;
    score(score>1) = 1;
end