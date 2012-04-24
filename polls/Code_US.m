% import the number of positives and negatives and count the sentiment
% score which is defined as the ratio of positive versus 
% negative messages on the topic on each day, the total number of positive
% words is 1389 and the number of negative words is 2162

% Consumer confidence considers Economy, Job, and Jobs
% Presidential approval considers Obama
% Elections considers Obama and Mccain

% The sentiment ratio for Economy
Data1 = importdata('C:\Data\Research\NLP\MATLAB\Code\JobsPositives.txt');
Data2 = importdata('C:\Data\Research\NLP\MATLAB\Code\JobsNegatives.txt');

Positives = Data1.data;
Negatives = Data2.data;
Ratio = Positives./Negatives;
% Ratio = Ratio(~isnan(Ratio));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth the sentiment ratio by computing the simple moving average over a 
% window of the past k days, where k indicates the number of previous 
% data points used with the current data point when calculating the moving average. 
k = 1; % no smoothing
AverageRatio_k = tsmovavg(Ratio, 's', k, 1) ;
% We use interpolation to deal with missing data
AverageRatio1 = interp1(1:161,AverageRatio_k,1:161,'cubic'); 
subplot(2,2,1),plot(1:161, AverageRatio_k, 1:161, AverageRatio1, 'r')
% subplot(2,2,1),plot(AverageRatio_k, 'r')
axis([0 80 0 2])
xlabel('No smoothing')
ylabel('Sentiment Ratio')
hold on

% Plot Figure 5 to exam the effects of using different windows on Jobs,
% i.e., k=3, k=7, k=15
k1 = 3; 
AverageRatio_k1 = tsmovavg(Ratio, 's', k1, 1) ;
AverageRatio2 = interp1(1:161,AverageRatio_k1,1:161,'cubic'); 
subplot(2,2,2),plot(1:161, AverageRatio_k1, 1:161, AverageRatio2, 'b')
% subplot(2,2,2),plot(AverageRatio1, 'b')
axis([0 80 0 2])
xlabel('Past 3 days')
ylabel('Sentiment Ratio')

k2 = 5;
AverageRatio_k2 = tsmovavg(Ratio, 's', k2, 1) ;
AverageRatio3 = interp1(1:161,AverageRatio_k2,1:161,'cubic'); 
subplot(2,2,3),plot(1:161, AverageRatio_k2, 1:161, AverageRatio3, 'black')
% subplot(2,2,3),plot(AverageRatio2, 'black')
axis([0 80 0 2])
xlabel('Past 5 days')
ylabel('Sentiment Ratio')

k3 = 7;
AverageRatio_k3 = tsmovavg(Ratio, 's', k3, 1) ;
AverageRatio4 = interp1(1:161,AverageRatio_k3,1:161,'cubic'); 
subplot(2,2,4),plot(1:161, AverageRatio_k3, 1:161, AverageRatio4, 'g')
% subplot(2,2,4),plot(AverageRatio3, 'g')
axis([0 80 0 2])
xlabel('Past week')
ylabel('Sentiment Ratio')
% xlabel('From 24/07/2009 to 31/12/2009')
% ylabel('Sentiment Ratio')
% legend('No smoothing', 'Past 3 days', 'Past 5 dys', 'Past week')
hold off

% Compare the trend using sentiment ratio with that of poll data based on
% the original data (k=1) and smoothed data (k=3)
% [num, txt, raw] = xlsread('DataSet-ECONINDEX122840-20120126')
[num, txt, raw] = xlsread('DataSet-Economy_161')
PollData = num(:,2);
% b = num(:,1);
n = length(PollData);
figure, 
subplot(2,1,1),plot(1:161, AverageRatio_k, 1:161, AverageRatio1, 'r')
hold on
plot(1:161, AverageRatio_k1, 1:161, AverageRatio2, 'b')
ylabel('Sentiment Ratio')
axis([0 80 0 2])
legend('k=3', 'k=1')
hold off
subplot(2,1,2), plot(1:80,PollData(1:80),'b')
xlabel('From 24/07/2009 to 11/10/2009')
ylabel('Gallup Economic Confidence Index')
% plot(1:n,b)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the correlation betweent the sentiment ratio and poll results by
% introduing a lag hyperparameter L (-30:30) to plot Figure 7.
% The L = k positions are marked on each curve
% Text-poll cross-correlation plots using a linear least-squares model
% where y is the poll outcome and x is the sentiment ratio
% a is the slope and b is the bias, which refer to the equations in the
% paper
[num, txt, raw] = xlsread('DataSet-Economy_161.xls');
cmap = colormap;
Tag = 1;
Day = [];
Ratio = [];
Poll = [];

for i = 1:length(Positives)
    if Positives(i)>0 & Negatives(i)>0
       Day(Tag) = i;
       Ratio(Tag) = Positives(i)./Negatives(i);
       Poll(Tag) = num(i,2);
       Tag = Tag+1;
    end
end

figure,
hold on

% We use the sentiment ratio with k from 1 to 16
for k = 1:2:16
    c = [];
    x = zeros(1,length(Ratio));
    for i = k:length(Ratio)
        x(i) = sum(Ratio(i-k+1:i))/k;
    end
    LTag = 1;
    for L = -30:30
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
            if x(i)>0 & y(i)<0
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
    plot(-30:30,c,'color',cmap(k*4,:));
    plot([-k -k],[c(30-k)-0.1 c(30-k)+0.1],'color',cmap(k*4,:),'linewidth',2);
    xlabel('Poll Lag with Different Smoothing Windows');
    ylabel('Correlation against Gallup')
end
