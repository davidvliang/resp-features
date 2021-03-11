% Uses JADE to find source signals from two sources in IQ modulated form
% uses IQ_ICA_func for computation

% Import data
load('data.mat')
I1 = data2file(:,1);
Q1 = data2file(:,2);
I2 = data2file(:,3);
Q2 = data2file(:,4);
CB1 = data2file(:,5);
CB2 = data2file(:,6);
time = data2file(:,7);

% find indices of set intervals
interval=10; % find intervals every 10 seconds
numIndices = 1;
indices = [];
for i=1:length(time)
    if (time(i) > numIndices*interval)
        indices(numIndices) = i;
        numIndices = numIndices + 1;
    end
end

% choose which interval to analyze
interval_choice = 10;
% calculate the range of values for that choice
if (interval_choice == 0)
    r = 1:indices(1);
else
    r = indices(interval_choice):indices(interval_choice+1);
end

% demodulate, ICA, and normalize
[source1, source2] = IQ_ICA_func(I1(r), Q1(r), I2(r), Q2(r));

% plot chestband
plotCB = 0;
if (plotCB)
    hold on;
    figure();
    plot(time(r), normalize(CB1(r)));
    plot(time(r), normalize(CB2(r)));
    legend('cb1', 'cb2');
    xlim([time(r(1)) time(r(end))]);
end

% show the results
plotSource = 1;
if (plotSource)
    figure();
    hold on;
    plot(time(r), -source1);
    plot(time(r), normalize(CB1(r)));
    plot(time(r), -source2);
    plot(time(r), normalize(CB2(r)));
    legend('src1', 'cb1', 'src2', 'cb2');
    xlim([time(r(1)) time(r(end))]);
end
