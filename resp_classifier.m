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
label_definition = {1 , "src1"; 2, "src2";};

% Gather All Features
X = [
        fwpt_feat_1, morph_feat_1;
        fwpt_feat_2, morph_feat_2
    ];

y = [
        ones(size(fwpt_feat_1,1),1); 
        2.*ones(size(fwpt_feat_2,1),1)
    ];

% Split Training and Test Data
[X_train, X_test, y_train, y_test] = SplitTrainTest(X, y, .70);


%% Support Vector Machine (SVM)

% model = fitcsvm(X_train, y_train, 'KernelFunction', 'polynomial', 'PolynomialOrder', 1);
model = fitcsvm(X_train, y_train, 'KernelFunction', 'linear');


%% Performance Evaluation and ROC Curve

p = predict(model, X_test);

fprintf('Test Accuracy: %f\n', mean(double(p == y_test)) * 100);
