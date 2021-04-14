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

%% Preprocess Data
% start = find(time==201, 1);
start = 1;
stop = find(time==500, 1); % define the stopping time (default: 200seconds)
% stop = length(time);

% splice data based on start/stop points
I1 = I1(start:stop);
Q1 = Q1(start:stop);
I2 = I2(start:stop);
Q2 = Q2(start:stop);
% CB1 = CB1(start:stop);
% CB2 = CB2(start:stop);
time = time(start:stop);


%% ICA

% demodulate, ICA, and normalize
[source1, source2] = IQ_ICA_func(I1, Q1, I2, Q2);

% Flip signals (this is done so it matches the chest band)
source1 = -source1;
source2 = -source2;

%% EMD
[src1_filtered, k_src1, imf_src1, ~] = filter_by_emd(source1);
[src2_filtered, k_src2, imf_src2, ~] = filter_by_emd(source2);

%% Segmentation
minpeakdist = 200;
minpeakprom = 0.5;

% Source 1
[max1, max1_loc] = findpeaks(src1_filtered, 'MinPeakDistance', minpeakdist, 'MinPeakProminence', minpeakprom);
[min1, min1_loc] = findpeaks(-src1_filtered, 'MinPeakDistance', minpeakdist, 'MinPeakProminence', minpeakprom);

% Source 2
[max2, max2_loc] = findpeaks(src2_filtered, 'MinPeakDistance', minpeakdist, 'MinPeakProminence', minpeakprom);
[min2, min2_loc] = findpeaks(-src2_filtered, 'MinPeakDistance', minpeakdist, 'MinPeakProminence', minpeakprom);

% Line up half-segments
peaks1 = clean_segments(max1, max1_loc, min1, min1_loc);
peaks2 = clean_segments(max2, max2_loc, min2, min2_loc);   

% Plot Peaks
figure(1);
subplot(2,1,1);
hold on;
plot(time(peaks1(:,2)), peaks1(:,1), 'linestyle','none', 'Marker', 'o', 'Color', 'b'); 
plot(time(peaks1(:,4)), peaks1(:,3), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
plot(time, src1_filtered, 'Color', 'k');
title('source1');
hold off;
subplot(2,1,2);
hold on;
plot(time(peaks2(:,2)), peaks2(:,1), 'linestyle','none', 'Marker', 'o', 'Color', 'b'); 
plot(time(peaks2(:,4)), peaks2(:,3), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
plot(time, src2_filtered, 'Color', 'k');
title('source2');
hold off;


%% Morphological Features

figure(2);
% pause;

% initialize feature vectors
morph_feat_1 = zeros(size(peaks1,1)-1, 100); % 100 feats/segment

% trough
Ht =  0; % height of the wave trough
Ht_p = 0; % p% of height H_t from minimum point
It_p = 0; % intercept of the waveform at height Ht_p

% crest
Hc =  0; % height of the wave crest
Hc_p = 0; % p% of height H_t from maximum point
Ic_p = 0; % intercept of the waveform at height Hc_p


f_idx = 1;
for i = 1:1:size(peaks1,1)-1 % for each segment (later consider overlappin segments)
    
    % Figure label
    clf;
    sgtitle(['Source 1: ', 'Segment ', num2str(f_idx), ' of ', num2str(size(morph_feat_1,1))]);
    
    % TROUGH
    xt = time(peaks1(i,2):peaks1(i+1,2),1);
    yt = src1_filtered(peaks1(i,2):peaks1(i+1,2),1);    
    Ht = min(peaks1(i,1), peaks1(i+1,1)) - peaks1(i,3); %modified to accomodate for uneven values
   
    t_idx = 1;
    for p = .02:.02:1
        
        % define height percentage (Ht_p)
        Ht_p = Ht * p + peaks1(i,3);
        if p == 1   % just to avoid sigfig issue
            Ht_p = min(peaks1(i,1), peaks1(i+1,1)); 
        end
        
        % find intersections
        [x0, y0] = intersections(xt, yt, xt, Ht_p*ones(size(xt)));
        x0_greater = min(x0(x0>time(peaks1(i,4))));
        x0_lesser = max(x0(x0<time(peaks1(i,4))));
        
        % Compute width
        It_p = x0_greater - x0_lesser;

        % visualize Ht_p
        figure(2);
        subplot(1,2,1);
        title('trough');
        hold on;
        plot(xt, yt, 'Color', 'black');
        plot(xt, ones(size(xt))*Ht_p, 'Color', 'm'); 
        plot(time(peaks1(i,4)), src1_filtered(peaks1(i,4)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([x0_lesser, x0_greater], Ht_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
    
        % morphological features
        morph_feat_1(f_idx, t_idx) = It_p / (Ht*p);
        
        % increment feature index
        t_idx = t_idx + 1;

    end
    
   
    % CREST
    xc = time(peaks1(i,4):peaks1(i+1,4),1);
    yc = src1_filtered(peaks1(i,4):peaks1(i+1,4),1);    
    Hc = peaks1(i+1,1) - max(peaks1(i,3), peaks1(i+1,3)); %modified to accomodate for uneven values
   
    c_idx = 51;
    for p = .02:.02:1
        
        % define crest percentage (Hc_p)
        Hc_p = peaks1(i+1,1) - Hc * p;
        if p == 1   % just to avoid sigfig issue
            Hc_p = max(peaks1(i,3), peaks1(i+1,3)); 
        end
        
        % find intersections
        [x0, y0] = intersections(xc, yc, xc, Hc_p*ones(size(xc)));
        x0_greater = min(x0(x0>time(peaks1(i+1,2))));
        x0_lesser = max(x0(x0<time(peaks1(i+1,2))));
        
        % Compute width
        Ic_p = x0_greater - x0_lesser;

        % visualize Hc_p
        figure(2);
        subplot(1,2,2);
        title('crest');
        hold on;
        plot(xc, yc, 'Color', 'black');
        plot(xc, ones(size(xc))*Hc_p, 'Color', 'm'); 
        plot(time(peaks1(i+1,2)), src1_filtered(peaks1(i+1,2)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([x0_lesser, x0_greater], Hc_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
    

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        % morphological features
        morph_feat_1(f_idx, c_idx) = Ic_p / (Hc*p);
        
        % increment feature index
        c_idx = c_idx + 1;

    end
    

    % increment vector index
    f_idx = f_idx + 1;
    
end


% initialize feature vectors
morph_feat_2 = zeros(size(peaks2,1)-1, 100); % 100 feats/segment

% trough
Ht =  0; % height of the wave trough
Ht_p = 0; % p% of height H_t from minimum point
It_p = 0; % intercept of the waveform at height Ht_p

% crest
Hc =  0; % height of the wave crest
Hc_p = 0; % p% of height H_t from maximum point
Ic_p = 0; % intercept of the waveform at height Hc_p


f_idx = 1;
for i = 1:1:size(peaks2,1)-1 % for each segment (later consider overlappin segments)
    
    % Figure label
    clf;
    sgtitle(['Source 2: ', 'Segment ', num2str(f_idx), ' of ', num2str(size(morph_feat_2,1))]);
    
    % TROUGH
    xt = time(peaks2(i,2):peaks2(i+1,2),1);
    yt = src2_filtered(peaks2(i,2):peaks2(i+1,2),1);    
    Ht = min(peaks2(i,1), peaks2(i+1,1)) - peaks2(i,3); %modified to accomodate for uneven values
   
    t_idx = 1;
    for p = .02:.02:1
        
        % define height percentage (Ht_p)
        Ht_p = Ht * p + peaks2(i,3);
        if p == 1   % just to avoid sigfig issue
            Ht_p = min(peaks2(i,1), peaks2(i+1,1)); 
        end
        
        % find intersections
        [x0, y0] = intersections(xt, yt, xt, Ht_p*ones(size(xt)));
        x0_greater = min(x0(x0>time(peaks2(i,4))));
        x0_lesser = max(x0(x0<time(peaks2(i,4))));
        
        % Compute width
        It_p = x0_greater - x0_lesser;

        % visualize Ht_p
        figure(2);
        subplot(1,2,1);
        title('trough');
        hold on;
        plot(xt, yt, 'Color', 'black');
        plot(xt, ones(size(xt))*Ht_p, 'Color', 'm'); 
        plot(time(peaks2(i,4)), src2_filtered(peaks2(i,4)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([x0_lesser, x0_greater], Ht_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
    

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        % morphological features
        morph_feat_2(f_idx, t_idx) = It_p / (Ht*p);
        
        % increment feature index
        t_idx = t_idx + 1;

    end
    
   
    % CREST
    xc = time(peaks2(i,4):peaks2(i+1,4),1);
    yc = src2_filtered(peaks2(i,4):peaks2(i+1,4),1);    
    Hc = peaks2(i+1,1) - max(peaks2(i,3), peaks2(i+1,3)); %modified to accomodate for uneven values
   
    c_idx = 51;
    for p = .02:.02:1
        
        % define crest percentage (Hc_p)
        Hc_p = peaks2(i+1,1) - Hc * p;
        if p == 1   % just to avoid sigfig issue
            Hc_p = max(peaks2(i,3), peaks2(i+1,3)); 
        end
        
        % find intersections
        [x0, y0] = intersections(xc, yc, xc, Hc_p*ones(size(xc)));
        x0_greater = min(x0(x0>time(peaks2(i+1,2))));
        x0_lesser = max(x0(x0<time(peaks2(i+1,2))));
        
        % Compute width
        Ic_p = x0_greater - x0_lesser;

        % visualize Hc_p
        figure(2);
        subplot(1,2,2);
        title('crest');
        hold on;
        plot(xc, yc, 'Color', 'black');
        plot(xc, ones(size(xc))*Hc_p, 'Color', 'm'); 
        plot(time(peaks2(i+1,2)), src2_filtered(peaks2(i+1,2)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([x0_lesser, x0_greater], Hc_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
    

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        % morphological features
        morph_feat_2(f_idx, c_idx) = Ic_p / (Hc*p);
        
        % increment feature index
        c_idx = c_idx + 1;

    end
    

    % increment vector index
    f_idx = f_idx + 1;
    
end



%% Fuzzy Wavelet Packet Transform Algorithm (FWPT) Features

% User-defined Parameters
winsize = 256; % window size
wininc = 2; % window increment
J = 7; % number of decomposition levels (default: 7)

% Initialize feature vectors
fwpt_feat_1 = zeros(size(peaks1,1)-1, 510); % 510 feats/segment
fwpt_feat_2 = zeros(size(peaks2,1)-1, 510);

f_idx = 1;
for i = 1:1:size(peaks1,1)-1
    seg = src1_filtered(peaks1(i,2):peaks1(i+1,4),1);    
    subspace = getmswpfeat(seg,winsize,wininc,J,'matlab');
    fwpt_feat_1(f_idx,:) = [skewness(subspace, 1, 1) std(subspace, 0, 1)];
    f_idx = f_idx + 1;
end

f_idx = 1;
for i = 1:1:size(peaks2,1)-1
    seg = src2_filtered(peaks2(i,2):peaks2(i+1,4),1);    
    subspace = getmswpfeat(seg,winsize,wininc,J,'matlab');
    fwpt_feat_2(f_idx,:) = [skewness(subspace, 1, 1) std(subspace, 0, 1)];
    f_idx = f_idx + 1;
end


%% Save Features

% save(['fwpt_1_ws', num2str(winsize), '_wi', num2str(wininc), '.mat'], 'fwpt_feat_1');
% save(['fwpt_2_ws', num2str(winsize), '_wi', num2str(wininc), '.mat'], 'fwpt_feat_2');
% save('morph_1.mat', 'morph_feat_1');
% save('morph_2.mat', 'morph_feat_2');


