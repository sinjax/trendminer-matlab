% import the number of positives and negatives and count the sentiment
% score which is defined as the ratio of positive versus 
% negative messages on the topic on each day, the total number of positive
% words is 1389 and the number of negative words is 2162

% The sentiment ratio for Conservative, Labour, and Liberal Democrats
Data1 = importdata('C:\Data\Research\NLP\MATLAB\Code\results_con(pre).txt');
Data2 = importdata('C:\Data\Research\NLP\MATLAB\Code\results_lab(pre).txt');
Data3 = importdata('C:\Data\Research\NLP\MATLAB\Code\results_lib(pre).txt');

Con = Data1.data;
Ratio_Con = Con(:,3);
Lab = Data2.data;
Ratio_Lab = Lab(:,3);
Lib = Data3.data
Ratio_Lib = Lib(:,3);
plot(1:123, Ratio_Con, 'r')
hold on
plot(1:123, Ratio_Lab, 'g')
plot(1:123, Ratio_Lib)
% Ratio = Ratio(~isnan(Ratio));
axis([0 130 0 3])
xlabel('Days from 2010-01-01 to 2010-05-06')
ylabel('Sentiment Ratio')
legend('Con', 'Lab', 'Lib')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth the sentiment ratio by computing the simple moving average over a 
% window of the past k days, where k indicates the number of previous 
% data points used with the current data point when calculating the moving average. 
k = 1; % no smoothing
AverageRatioCon_k = tsmovavg(Ratio_Con, 's', k, 1) ;
AverageRatioLab_k = tsmovavg(Ratio_Lab, 's', k, 1) ;
AverageRatioLib_k = tsmovavg(Ratio_Lib, 's', k, 1) ;
% We use interpolation to deal with missing data
% AverageRatio1 = interp1(1:161,AverageRatio_k,1:161,'cubic'); 
figure
subplot(2,2,1),plot(1:123, AverageRatioCon_k, 'r')
hold on
plot(1:123, AverageRatioLab_k, 'g')
plot(1:123, AverageRatioLib_k, 'b')
% subplot(2,2,1),plot(AverageRatio_k, 'r')
axis([0 130 1 3])
xlabel('No Smoothing')
ylabel('Sentiment Ratio')
hold on

% Plot Figure 5 to exam the effects of using different windows on Jobs,
% i.e., k=3, k=7, k=15
k1 = 5; 
AverageRatioCon_k1 = tsmovavg(Ratio_Con, 's', k1, 1) ;
AverageRatioLab_k1 = tsmovavg(Ratio_Lab, 's', k1, 1) ;
AverageRatioLib_k1 = tsmovavg(Ratio_Lib, 's', k1, 1) ; 
subplot(2,2,2),plot(1:123, AverageRatioCon_k1, 'r')
hold on
plot(1:123, AverageRatioLab_k1, 'g')
plot(1:123, AverageRatioLib_k1, 'b')
% subplot(2,2,2),plot(AverageRatio1, 'b')
axis([0 130 1 3])
xlabel('Past 5 Days')
ylabel('Sentiment Ratio')

k2 = 7;
AverageRatioCon_k2 = tsmovavg(Ratio_Con, 's', k2, 1) ;
AverageRatioLab_k2 = tsmovavg(Ratio_Lab, 's', k2, 1) ;
AverageRatioLib_k2 = tsmovavg(Ratio_Lib, 's', k2, 1) ; 
subplot(2,2,3),plot(1:123, AverageRatioCon_k2, 'r')
hold on
plot(1:123, AverageRatioLab_k2, 'g')
plot(1:123, AverageRatioLib_k2, 'b')
% subplot(2,2,3),plot(AverageRatio2, 'black')
axis([0 130 1 3])
xlabel('Past 7 Days')
ylabel('Sentiment Ratio')

k3 = 14;
AverageRatioCon_k3 = tsmovavg(Ratio_Con, 's', k3, 1) ;
AverageRatioLab_k3 = tsmovavg(Ratio_Lab, 's', k3, 1) ;
AverageRatioLib_k3 = tsmovavg(Ratio_Lib, 's', k3, 1) ;
subplot(2,2,4),plot(1:123, AverageRatioCon_k3, 'r')
hold on
plot(1:123, AverageRatioLab_k3, 'g')
plot(1:123, AverageRatioLib_k3, 'b')
% subplot(2,2,4),plot(AverageRatio3, 'g')
axis([0 130 1 3])
xlabel('Past 14 Days')
ylabel('Sentiment Ratio')
hold off

% Compare the trend using sentiment ratio with that of poll data based on
% the original data (k=1) and smoothed data (k=3)
% [num, txt, raw] = xlsread('DataSet-ECONINDEX122840-20120126')
PollData = importdata('C:\Data\Research\NLP\MATLAB\Code\UK_Interpolation(pre).txt');
Poll_Con = PollData.data(:,1);
Poll_Lab = PollData.data(:,2);
Poll_Lib = PollData.data(:,3);
% b = num(:,1);
figure, 
subplot(2,1,1),plot(1:123, AverageRatioCon_k2, 'r')
hold on
plot(1:123, AverageRatioLab_k2, 'g')
plot(1:123, AverageRatioLib_k2, 'b')
ylabel('Sentiment Ratio')
axis([0 130 1 3])
legend('Con', 'Lab', 'Lib')
hold off
subplot(2,1,2), plot(1:123,Poll_Con,'r')
hold on
plot(1:123,Poll_Lab,'g')
plot(1:123,Poll_Lib,'b')
xlabel('Days from 2010-01-04 to 2010-05-06')
ylabel('YouGov Polling Data (Percentage)')
legend('Con', 'Lab', 'Lib')
% plot(1:n,b)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the correlation betweent the sentiment ratio and poll results by
% introduing a lag hyperparameter L (-20:20) to plot Figure 7.
% The L = k positions are marked on each curve
% Text-poll cross-correlation plots using a linear least-squares model
% where y is the poll outcome and x is the sentiment ratio
% a is the slope and b is the bias, which refer to the equations in the
% paper
% [num, txt, raw] = xlsread('DataSet-Economy_161.xls');

% For Conservative party
cmap = colormap;
Tag = 1;
Day = [];
Ratio = [];
Poll = [];

figure,
hold on

for i = 1:length(Poll_Con)
    if Poll_Con(i)>0
       Day(Tag) = i;
       Ratio(Tag) = Ratio_Con(i);
       Poll(Tag) = Poll_Con(i);
       Tag = Tag+1;
    end
end

% We use the sentiment ratio with k from 1 to 16
for k = 1:2:16
    c = [];
    x = zeros(1,length(Ratio));
    for i = k:length(Ratio)
        x(i) = sum(Ratio(i-k+1:i))/k;
    end
    LTag = 1;
    for L = -20:20
        y = zeros(1,length(Poll));
        if L<0
            le = length(Poll)+L;
            for i = 1:le
                y(i-L) = Poll(i);
            end
        else
            le = length(Poll)-L;
            for i = 1:le
                y(i) = Poll(i+L);
            end
        end
        Tag = 1;
        xx = [];
        yy = [];
        for i = 1:length(Ratio)
            if x(i)>0 & y(i)>0
                xx(Tag) = x(i);
                yy(Tag) = y(i);
                Tag = Tag+1;
            end
        end
%         figure
%         subplot(1,2,1);
%         plot(xx,yy,'o');
%         subplot(1,2,2);
%         hold on
%         plot(xx/mean(xx),'r');
%         plot(yy/mean(yy),'b');
%         hold off
        
        [r p] = corrcoef(xx,yy);
        
        c(LTag) = r(1,2);
        LTag = LTag+1;        
    end
    plot(-20:20,c,'color',cmap(k*4,:));
    plot([-k -k],[c(30-k)-0.1 c(30-k)+0.1],'color',cmap(k*4,:),'linewidth',2);
    xlabel('Poll Lag with Different Smoothing Windows for Con');
    ylabel('Correlation against YouGov')
end

