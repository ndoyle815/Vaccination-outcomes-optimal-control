% script to run simulations of control model with and without vaccination,
% produces Figure 1 (cumulative vaccinations) and Figure 3 (simulations) from the report
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

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 800;  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 100 150 600; 50 100 150 200; 100 300 500 600; 250 350 425 500];
strategies = [1:length(thresholds)];

% plotting preperation
figure('Position',[200 400 900 900])
stIDs = {'S1', 'S2', 'S3', 'S4'};
stnames = {'(Cautious easing)', '(Suppression)', '(Slow control)', '(Rapid control)'};
stpos = [50 200 350 600; 50 200 350 500; 100 300 475 625; 175 325 475 625];

tic
for strat = strategies
    % Preliminary run - no control, 30 day build-up
    para0.vstart = vstarts;
    [Prelim, Prelim_ICs] = Get_ICs(para0);
        
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
    [Classes, burden, stringency, peak_hospital] = ODEmodel(para,Prelim_ICs); 

    % Post-Process for epidemic metrics
    [Classes, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, nx, ix1, ix2, ~, ~] = PostProcessor(Classes);

    plotidx = 2*strat - 1;  % index for subfigure

    % plotting active hospitalisations
    subplot(length(strategies),2,plotidx)

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

    if plotidx > 2*(length(strategies)-1)
        xlabel('Time (days)')
    end

    ylabel({char(stIDs(strat)) ; char(stnames(strat))})

    if plotidx == 1
        title('Active $I^H(t)$');
    end

    axis([0 maxtime 0 1300])

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

    axis([0 maxtime 0 1300])
        
    grid on
 

    % plotting cumulative hospitalisations
    plotidx = plotidx + 1;
    subplot(length(strategies),2,plotidx)

    for i = ix1'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 40000 40000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    for i = ix2'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 40000 40000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.Hosp,2)/1000, 'k', 'LineWidth', 2.5)
    
    if vstarts < maxtime
        if plotidx == 2
            xline(para.vstart,'-','Vaccine Arrival','Color','b','Linewidth',2,'FontSize',16,'Interpreter','latex','LabelOrientation','horizontal')
        else
            xline(para.vstart,'-','Color','b','Linewidth',2)
        end
    end
        
    axis([0 maxtime 0 max([1.01*sum(Classes.Hosp(end,:))/1000,30])])
    
    if plotidx > 2*(length(strategies)-1)
        xlabel('Time (days)')
    end

    if plotidx == 2
        title('Cumulative $I^H(t)$ (thousands)');
    end

    set(gca,'FontSize',16)
    grid on

end
toc

%save figure
saveas(gcf,strcat('./sim_images/','simulation_',num2str(para.vstart),'.png'))


% Plotting cumulative vaccinations
f = figure(2);
f.Position = [1250 400 450 300];
plot(Classes.t, Classes.V(:,1)./1000, 'LineWidth', 2.5, 'DisplayName', '0-19')
hold on
plot(Classes.t, Classes.V(:,2)./1000, 'LineWidth', 2.5, 'DisplayName', '20-64')
hold on
plot(Classes.t, Classes.V(:,3)./1000, 'LineWidth', 2.5, 'DisplayName', '65+')
hold on
plot(Classes.t, sum(Classes.V,2)./1000, 'k', 'LineWidth', 2.5, 'DisplayName', 'Total')

axis([0 maxtime 0 sum(para.N)./1000])
xline(para.vstart,'--','DisplayName','Arrival','Color','k')
xlabel('Time (days)')
ylabel('Population (thousands)')
title('Cumulative Vaccinations')
legend('Interpreter','latex','Location','west')
grid on

saveas(gcf,'./vacc_images/Tvacc.png')
