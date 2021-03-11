%% Clear Workspace
clc;
close all;
clear all;

%% Import Data
load('Data/data.mat')
I1 = data2file(:,1);
Q1 = data2file(:,2);
I2 = data2file(:,3);
Q2 = data2file(:,4);
% CB1 = data2file(:,5);
% CB2 = data2file(:,6);
time = data2file(:,7);

%% ICA
stop = find(time==100, 1); % define the stopping time (default: 100seconds)

% demodulate, ICA, and normalize
[source1, source2] = IQ_ICA_func(I1(1:stop), Q1(1:stop), I2(1:stop), Q2(1:stop));


%% Segmentation
[max1, max_loc1] = findpeaks(source1(1:stop),'MinPeakDistance', 20);
[min_1, min_loc1] = findpeaks(-source1(1:stop),'MinPeakDistance', 20);

figure;
hold on;
% plot(max_loc1, max1, 'o'); % doesn't work
plot(time(1:stop), -source1(1:stop));
plot(time(1:stop), -source2(1:stop));
legend('src1', 'src2');
hold off;




%% Morphological Features




%% Fuzzy Wavelet Packet Transform Algorithm (FWPT) Features

% Features = getmswpfeat(x,winsize,wininc,J,toolbox); %unfinished


%% Format features (70-30 split)

% Define Labels
label_definition = {1 , "src1"; 2, "src2";};

% Split Training and Test Data
[X_train, X_test, y_train, y_test] = SplitTrainTest(X, y, .70);


% Support Vector Machine (SVM)



%% Performance Evaluation and ROC Curve


