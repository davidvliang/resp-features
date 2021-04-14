%% Clear Workspace 

% close all;
clear all;
clc;


%% Load Features

load('fwpt_1_ws256_wi2.mat');
load('fwpt_2_ws256_wi2.mat');
load('morph_1.mat');
load('morph_2.mat');


%% Format features (70-30 split) 

% Define Labels
labels = ["src1"; "src2"];

% Gather All Features
X = [
        fwpt_feat_1, morph_feat_1;
        fwpt_feat_2, morph_feat_2
    ];

y = [
        ones(size(fwpt_feat_1,1),1); 
        2.*ones(size(fwpt_feat_2,1),1)
    ];

%% Support Vector Machine (SVM)
runs = 100;

test_accuracies = zeros(runs,1);

for i = 1:runs

    % Split Training and Test Data
    [X_train, X_test, y_train, y_test] = SplitTrainTest(X, y, .70);


    % Train Model
    model = fitcsvm(X_train, y_train, 'KernelFunction', 'polynomial', 'PolynomialOrder', 2);
%     model = fitcsvm(X_train, y_train, 'KernelFunction', 'linear');

    % Performance Evaluation and ROC Curve
    [p, score] = predict(model, X_test);
%     table(labels(y_test),labels(p), max(score,[],2),'VariableNames', {'TrueLabels','PredictedLabels','Score'});
    test_accuracies(i) = mean(double(p == y_test)) * 100;
    fprintf('Test Accuracy: %f\n', test_accuracies(i));

end

fprintf('\n--\nMean Test Accuracy: %f\n', mean(test_accuracies));
