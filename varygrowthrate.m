% scenario where control strategies are simulated with variable recovery
% rate gamma, but R0 is rescaled so it is equal to 3.0 in all cases
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination
vstarts = 2160;

% gammas and R0 factors
gammas = [0.4641, 1/5.6, 1/12];
betamat = [1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];
R0scales = [2.443, 0.945, 0.442];
E0s = [0.45e-4, 2e-4, 8e-4];
ngams = length(gammas);

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 800;  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 100 150 600; 50 100 150 200; 100 300 500 600; 250 350 425 500];
strategies = [1:length(thresholds)];

% plotting preperation
figure('Position',[200 400 500*ngams 1000])
stIDs = {'S1', 'S2', 'S3', 'S4'};
stnames = {'(Cautious easing)', '(Suppression)', '(Slow control)', '(Rapid control)'};
stpos = [50 200 350 600; 50 200 350 500; 100 300 475 625; 175 325 475 625];

tic
for gs = 1:ngams
    para0.gamma = gammas(gs);
    para0.zeta = gammas(gs);
    para0.beta = R0scales(gs).*betamat;
    para0.E0 = E0s(gs);

    Get_R0(para0)

for strat = strategies
    % Preliminary run - no control, 30 day build-up
    para0.vstart = vstarts;
    [Prelim, Prelim_ICs] = Get_ICs(para0);

    sum(Prelim_ICs.Hosp(end,:))
        
    % add control thresholds defined by strategy
    para = para0;
    para.maxtime = maxtime;
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.T12
        para.init = 1;
    else
        para.init = 2;
    end

    % Run model
    [Classes] = ODEmodel(para,Prelim_ICs); 

    % Post-Process for epidemic metrics
    [Classes, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, nx, ix1, ix2, ~, ~] = PostProcessor(Classes);

    plotidx = ngams*(strat-1) + gs;  % index for subfigure

    % plotting active hospitalisations
    subplot(length(strategies),ngams,plotidx)

    yyaxis left
    ax1 = gca;
    ax1.YColor = 'k';
    ax1.FontSize = 16;
    ax1.FontSizeMode = 'manual';

    for i = ix1'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    for i = ix2'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.IH,2), 'k', 'LineWidth', 2.5)
    
    if vstarts < maxtime
        if plotidx == 1
            xline(para.vstart,'-','Vaccine Arrival','Color','b','Linewidth',2,'FontSize',16,'Interpreter','latex','LabelOrientation','horizontal')
        else
            xline(para.vstart,'-','Color','b','Linewidth',2)
        end
    end

    if plotidx > ngams*(length(strategies)-1)
        xlabel('Time (days)')
    end

    if gs == 1
        ylabel({char(stIDs(strat)) ; char(stnames(strat))})
    end

    if strat == 1
        title(strcat('$\gamma = ',num2str(gammas(gs)),'$'));
    end

    axis([0 maxtime 0 1400])

    yyaxis right
    ax2 = gca;
    ax2.YColor = 'k';
    ax2.TickDir = 'none';
    ax2Y = ax2.YAxis(2,1);
    ax2Y.FontSize = 14;
    ax2.FontSizeMode = 'manual';

    plot(Classes.t, para.T01.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T12.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T10.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T21.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    yticks(stpos(strat,:))
    
    if para.T21 > para.T01
        yticklabels({'$T_{10}$','$T_{01}$','$T_{21}$','$T_{12}$'})
    else
        yticklabels({'$T_{10}$','$T_{21}$','$T_{01}$','$T_{12}$'})
    end

    axis([0 maxtime 0 1400])
        
    grid on
 

    

end
end
toc

%save figure
saveas(gcf,'./sim_images/simulation_vargrowthrate.png')

