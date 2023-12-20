clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')

% load cost outputs
load('./mats/cost_tensor.mat');

% load joint distributions
load('./mats/jointdists.mat');

% strategy markers
markers = {'o','^','s','d'};
cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

%% Computing expectations and plotting

strategies = [1:size(fs,2)];
weights = [0.0:0.01:1.0];
etas = [0.0:0.05:1.0];
Ts = [60:60:1080];
dweights = [1:5:size(fs,1)];  % only plot some weights with markers

% number of distributions
ndists = size(P,1);

% upper and lower percentiles for credible intervals
credintsize = 0.95;
LB = (1 - credintsize)/2;
UB = 1 - LB;

% find min cost (=1) or max cost (=2)
minormax = 1;

for dist = 1:ndists

    % matrix to store expected cost
    Ecosts = zeros(length(weights),length(strategies));

    % ..and lower / upper bounds of credible interval
    CostLB = 0.2.*ones(length(dweights),length(strategies));
    CostUB = 1.0.*ones(length(dweights),length(strategies));

    f = figure(dist);
    if dist == 1
        hht = 287.5;
    else
        hht = 250;
    end
    f.Position = [400 400 1250 hht];
    
    % obtain joint distribution
    PD = reshape(P(dist,:,:), [length(Ts) length(etas)]);

    % matrix to store probabilities each strategy is optimal per weight
    Poptimal = zeros(length(weights),length(strategies));

    dc = 0.0001;
    dws = [1, 1, 1, 1];

    for w = 1:length(weights)

        if minormax == 1
            [~,ix] = min([fs(w,1,:,:);fs(w,2,:,:);fs(w,3,:,:);fs(w,4,:,:)]);
        else
            [~,ix] = max([fs(w,1,:,:);fs(w,2,:,:);fs(w,3,:,:);fs(w,4,:,:)]);
        end
        ix = reshape(ix, [length(Ts) length(etas)]);

        for strat = strategies
            
            % obtain function values
            FWS = reshape(fs(w,strat,:,:), [length(Ts) length(etas)]);

            % expected cost
            Ecosts(w,strat) = sum(PD.*FWS,'all');

            opt = ix == strat;
            Poptimal(w,strat) = sum(PD.*opt,'all');

            if ismember(w,dweights)

                while sum(PD(FWS < CostUB(dws(strat),strat))) > UB
                    CostUB(dws(strat),strat) = CostUB(dws(strat),strat) - dc;
                end
                while sum(PD(FWS < CostLB(dws(strat),strat))) < LB
                    CostLB(dws(strat),strat) = CostLB(dws(strat),strat) + dc;
                end
                dws(strat) = dws(strat) + 1;
            end

        end
    end

    subplot(1,4,1)

    ax = gca;
    ax.InnerPosition(1) = ax.InnerPosition(1) - 0.05;
    if dist == 1
        ax.InnerPosition(2) = 0.2507*250/hht;
        ax.InnerPosition(4) = 0.6671*250/hht;
    else
        ax.InnerPosition(2) = 0.2507;
        ax.InnerPosition(4) = 0.6671;  
    end

    imagesc(etas, Ts, PD)
    xlab = xlabel('$\eta$');
    xlab.Position(2) = 1300;
    ylab = ylabel('$T$','Rotation',0);
    ylab.Position(1) = ylab.Position(1) - 0.06;
    xticks(etas(1:5:end))
    xtickangle(0)
    yticks(Ts(2:3:end))
    set(gca,'FontSize',18)

    subplot(1,4,[2 3])

    for strat = strategies
        plot(weights(dweights),Ecosts(dweights,strat),"Color",'k','LineWidth',1.5);
        hold on
        plot(weights(dweights(2:end-1)),Ecosts(dweights(2:end-1),strat),"Color",'k','Marker',markers{strat},'MarkerSize',8, ...
             'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5);
        hold on
        patch([weights(dweights) fliplr(weights(dweights))], [(CostLB(:,strat))' fliplr((CostUB(:,strat))')], ...
              cols(strat,:),'EdgeColor','none','FaceAlpha',0.4)
        hold on
    end

    axis([min(weights)+0.0 max(weights)-0.0 0.2 1])

    ax = gca;
    if dist == 1
        ax.InnerPosition(2) = 0.2507*250/hht;
        ax.InnerPosition(4) = 0.6671*250/hht;
    else
        ax.InnerPosition(2) = 0.2507;
        ax.InnerPosition(4) = 0.6671;        
    end

    ylabel('Cost')
    xlab = xlabel('$w$');
    xlab.Position(2) = 0.05;
    set(gca,'TickLength',[0 0])

    xticks([min(weights):0.25:max(weights)])
    yticks([0:0.2:1])
    set(gca,'FontSize',18)
    grid on

    if dist == 1
        leg = legend({'','S1 (Cautious Easing)','','','S2 (Suppression)','','','S3 (Slow Control)','','','S4 (Rapid Control)',''},'Interpreter','Latex',...
              'FontSize',18,'Orientation','horizontal','Position',[0.45 0.9 0.1 0.06]);
    end

    subplot(1,4,4)

    for strat = strategies
        plot(weights,Poptimal(:,strat),"Color",cols(strat,:),"LineWidth",2.5)
        hold on
    end

    axis([min(weights) max(weights) 0 1])

    ax = gca;
    ax.InnerPosition(1) = ax.InnerPosition(1) + 0.05;
    if dist == 1
        ax.InnerPosition(2) = 0.2507*250/hht;
        ax.InnerPosition(4) = 0.6671*250/hht;
    else
        ax.InnerPosition(2) = 0.2507;
        ax.InnerPosition(4) = 0.6671;  
    end

    ylabel('Probability')
    xlab = xlabel('$w$');
    xlab.Position(2) = -0.18;
    xticks([min(weights):0.25:max(weights)])
    xtickangle(0)

    set(gca,'FontSize',18)
    set(gca,'TickLength',[0 0])
    grid on

    saveas(f,strcat('./vacc_images/expcost_jointdists',num2str(dist),'.png'))

    %close all

end

