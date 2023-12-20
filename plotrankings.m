% script to merge subplots for control strategy rankings without
% vaccination into Figure 4 from report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',22)

A1 = load('./mats/StratRankings_1500_500.mat');
A2 = load('./mats/StratRankings_1500_1000.mat');
A3 = load('./mats/StratRankings_1500_1500.mat');
A4 = load('./mats/StratRankings_1075_1000.mat');

Ranks(:,:,1) = A1.rank;
Ranks(:,:,2) = A2.rank;
Ranks(:,:,3) = A3.rank;
Ranks(:,:,4) = A4.rank;
Hc = [1500, 1500, 1500, 1075];
maxtime = [500, 1000, 1500, 1000];

weights = [0:0.01:1];

cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

figure('Position',[400 400 1600 600])

for k = 1:4
    subplot(2,2,k)
    hold on
    for strat = 1:4
        plot(weights,Ranks(strat,:,k),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
    end
    hold off
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    yticks([1,2,3,4])
    ylabel('Ranking')
    
    if k >= 3
        xlabel('$w_1 \; (w_2 = 1 - w_1)$')
    end
    if k == 5
        legend({'S1 (Cautious Easing)','S2 (Suppression)','S3 (Slow Control)','S4 (Rapid Control)'},'Fontsize',14,'Interpreter','Latex','Location','east');
    end
    title(strcat('$H_c$ = ',{' '},num2str(Hc(k)),', ',' maxtime =  ',{' '},num2str(maxtime(k))))
    grid on
end

saveas(gcf,'./sim_images/Costoutputs_vacc.png')