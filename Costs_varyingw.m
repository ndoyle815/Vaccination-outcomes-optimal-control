% script to run control strategies without vaccination and rank according
% to the objective function, producing subplots for Figure 4 in the report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')

% Plotting
markers = {'o','^','s','d'};
cols = [0.9290 0.6940 0.1250; 0.3290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

% load default parameters
para0 = load('./mats/Parameters.mat');

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 1080;  % main simulation
vstarts = [2160, 360]; % vaccine times

% define strategy numbers and switching thresholds
thresholds = [50 100 150 600; 50 100 150 200; 100 300 500 600; 250 350 425 500];
strategies = [1:length(thresholds)];

% add control thresholds defined by strategy
para = para0;
para.init = 1;
para.maxtime = maxtime;
para.Hmax = 1499;        % modify hospital capacity

% define functional weights
weights = [0:0.01:1];  % varying w
w2 = 2;

f1 = figure();
if para.maxtime < 500
    hht = 250;
else
    hht = 200;
end
f1.Position = [600 600 600*length(vstarts) hht];

tic
for vs = 1:length(vstarts)
    para.vstart = vstarts(vs);
    
    % stores cost function outputs
    ns = length(strategies);
    nw = length(weights);
    fs = zeros(ns,nw);

    for strat = strategies
        % set switching thresholds
        para.T10 = thresholds(strat,1);
        para.T01 = thresholds(strat,2);
        para.T21 = thresholds(strat,3);
        para.T12 = thresholds(strat,4);

        % run preliminary simulation to get ICs
        [Prelim, Prelim_ICs] = Get_ICs(para0);

        % Run main simulation
        [Classes, burden, stringency, peak_hospital] = ODEmodel(para, Prelim_ICs);
    
        for w = 1:nw
            % evaluate cost function
            fs(strat,w) = CostFunction([weights(w), w2], para, burden, stringency, peak_hospital, 0);
        end
    end

    % rank strategies
    [~,idx] = min(fs);
    [sorted, ii] = sort(fs);
    [~,rank] = sort(ii);

    subplot(1,2,vs)

    for strat = strategies
        plot(weights,rank(strat,:),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
        hold on
    end

    ax = gca;
    if para.maxtime > 500
        ax.Position(2) = ax.Position(2) + 0.25;
        ax.Position(4) = ax.Position(4) - 0.25;
    else
        ax.Position(2) = ax.Position(2) + 0.18;
        ax.Position(4) = 0.603*200/hht;
    end

    set(gca, 'YDir','reverse')
    set(gca, 'FontSize',20)

    axis([min(weights) max(weights) 0.5 4.5])
    xticks(weights(1:20:end))
    yticks([1,2,3,4])
    xtickangle(0)
    xlab = xlabel('Weight $w$');

    if vs == 1
        ylab = ylabel('Ranking','Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
        ylab.Position(1) = ylab.Position(1) - 0.02;
    end
    grid on

end
toc

%if para.Hmax < 1200 && para.vstart < para.maxtime
if para.maxtime < 500
    legend({'S1 (Cautious Easing)','S2 (Suppression)','S3 (Slow Control)','S4 (Rapid Control)'},...
            'Interpreter','Latex','FontSize',18,'Orientation','horizontal','Position',[0.09 0.84 0.82 0.1])
end


saveas(gcf,strcat('./sim_images/CostFunction_',num2str(para.Hmax),'_',num2str(maxtime),'.png'))

% save(strcat('./mats/StratRankings_',num2str(para.vstart),'_',num2str(para.Hmax),'_',num2str(maxtime),'.mat'),"rank")
