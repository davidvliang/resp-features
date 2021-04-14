function [morph_matrix] = morph_features(time, input_sig, peak_vals)
%MORPH_FEATURES Given time, input signal, and peak values. Compute the morphological 
% features as described in [1]

figure(2);
% pause;

% initialize feature vectors
morph_matrix = zeros(size(peak_vals,1)-1, 100); % 100 feats/segment

% trough
Ht =  0; % height of the wave trough
Ht_p = 0; % p% of height H_t from minimum point
It_p = 0; % intercept of the waveform at height Ht_p

% crest
Hc =  0; % height of the wave crest
Hc_p = 0; % p% of height H_t from maximum point
Ic_p = 0; % intercept of the waveform at height Hc_p


f_idx = 1;
for i = 1:1:size(peak_vals,1)-1 % for each segment (later consider overlappin segments)
    
    % Figure label
    clf;
    sgtitle(['Source 1: ', 'Segment ', num2str(f_idx), ' of ', num2str(size(morph_matrix,1))]);
    
    % TROUGH
    xt = time(peak_vals(i,2):peak_vals(i+1,2),1);
    yt = input_sig(peak_vals(i,2):peak_vals(i+1,2),1);    
    Ht = min(peak_vals(i,1), peak_vals(i+1,1)) - peak_vals(i,3); %modified to accomodate for uneven values
   
    % CREST
    xc = time(peak_vals(i,4):peak_vals(i+1,4),1);
    yc = input_sig(peak_vals(i,4):peak_vals(i+1,4),1);    
    Hc = peak_vals(i+1,1) - max(peak_vals(i,3), peak_vals(i+1,3)); %modified to accomodate for uneven values
   
    p_idx = 1;
    for p = .02:.02:1
        
        % define trough height percentage (Ht_p)
        Ht_p = Ht * p + peak_vals(i,3);
        if p == 1   % just to avoid sigfig issue
            Ht_p = min(peak_vals(i,1), peak_vals(i+1,1)); 
        end
        
        % define crest height percentage (Hc_p)
        Hc_p = peak_vals(i+1,1) - Hc * p;
        if p == 1   % just to avoid sigfig issue
            Hc_p = max(peak_vals(i,3), peak_vals(i+1,3)); 
        end
        
        % find trough intersections
        [xint_t, ~] = intersections(xt, yt, xt, Ht_p*ones(size(xt)));
        xint_t_greater = min(xint_t(xint_t>time(peak_vals(i,4))));
        xint_t_lesser = max(xint_t(xint_t<time(peak_vals(i,4))));
        
        % find crest intersections
        [xint_c, ~] = intersections(xc, yc, xc, Hc_p*ones(size(xc)));
        xint_c_greater = min(xint_c(xint_c>time(peak_vals(i+1,2))));
        xint_c_lesser = max(xint_c(xint_c<time(peak_vals(i+1,2))));

        
        % Compute widths
        It_p = xint_t_greater - xint_t_lesser;
        Ic_p = xint_c_greater - xint_c_lesser;


        % visualize Ht_p
        figure(2);
        subplot(1,2,1);
        title('trough');
        hold on;
        plot(xt, yt, 'Color', 'black');
        plot(xt, ones(size(xt))*Ht_p, 'Color', 'm'); 
        plot(time(peak_vals(i,4)), input_sig(peak_vals(i,4)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([xint_t_lesser, xint_t_greater], Ht_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
        
        % visualize Hc_p
        figure(2);
        subplot(1,2,2);
        title('crest');
        hold on;
        plot(xc, yc, 'Color', 'black');
        plot(xc, ones(size(xc))*Hc_p, 'Color', 'm'); 
        plot(time(peak_vals(i+1,2)), input_sig(peak_vals(i+1,2)), 'linestyle','none', 'Marker', 'o', 'Color', 'r'); 
        plot([xint_c_lesser, xint_c_greater], Hc_p*ones(1,2), 'pg', 'MarkerSize',10);
        hold off;
        
        
        % morphological features
        morph_matrix(f_idx, p_idx) = It_p / (Ht*p);
        morph_matrix(f_idx, p_idx+50) = Ic_p / (Hc*p);

        
        % increment feature index
        p_idx = p_idx + 1;

    end
    
    % increment vector index
    f_idx = f_idx + 1;
    
end

end

