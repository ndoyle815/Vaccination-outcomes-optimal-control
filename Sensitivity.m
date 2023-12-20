% script to run simulations of control model with and without vaccination,
% produces Figure 1 (cumulative vaccinations) and Figure 3 (simulations) from the report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

grayColor = [.6 .6 .6];

% load default parameters
para0 = load('./mats/Parameters.mat');

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 800;  % main simulation

% Sensitivity parameter and selection
vars = {'\gamma','\delta','\rho','\tau','\omega','d_a','h_a'};
labs = {'gam','del','rho','tau','ome','da','ha'};

var_to_test = 2;

gammas = [0.1:0.02:0.3];
gammas(end+1) = 1/5.6;
rhos = [0.0:0.02:0.2];
rhos(end+1) = 0.1;
deltas = [0.05:0.02:0.25];
deltas(end+1) = 1/8.78;
taus = [0.1:0.03:0.4];
taus(end+1) = 0.25;
omegas = [0.0005:0.0002:0.0025];
omegas(end+1) = 1/750;
das = [0.5:0.1:1.5];
das(end+1) = 1;
has = [0.5:0.15:2.0];
has(end+1) = 1;

if var_to_test == 1
    var_arr = gammas;
elseif var_to_test == 2
    var_arr = deltas;
elseif var_to_test == 3
    var_arr = rhos;
elseif var_to_test == 4
    var_arr = taus;
elseif var_to_test == 5
    var_arr = omegas;
elseif var_to_test == 6
    var_arr = das;
elseif var_to_test == 7
    var_arr = has;
end

% define strategy numbers and switching thresholds
thresholds = [50 100 150 600; 50 100 150 200; 100 300 500 600; 250 350 425 500];
strategies = [1:length(thresholds)];

% store disease burden and stringency cost contributions
Burden = zeros(length(var_arr),length(strategies));
Stringency = zeros(length(var_arr),length(strategies));

% plotting preperation
figure('Position',[200 400 1000 1000])
stratnames = {'S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)'};
stpos = [50 200 350 600; 50 200 350 500; 100 300 475 625; 175 325 475 625];
sgtitle(strcat('Sensitivity, varying $',vars(var_to_test),' \in [',num2str(min(var_arr)),',',num2str(max(var_arr)),']$'),'FontSize',20)

tic
for idx = 1:length(var_arr)

    for strat = strategies
        % change sensitivity parameter
        if var_to_test == 1
            para0.gamma = gammas(idx);
            para0.zeta  = gammas(idx);
        elseif var_to_test == 2
            para0.delta = deltas(idx);
        elseif var_to_test == 3
            para0.rho = rhos(idx);
        elseif var_to_test == 4
            para0.tau = taus(idx);
        elseif var_to_test == 5
            para0.omega = omegas(idx);
        elseif var_to_test == 6
            para0.da = das(idx).*[0.0275; 0.1350; 0.5374];
        elseif var_to_test == 7
            para0.ha = has(idx).*[0.1268; 0.1173; 0.2430];
        end

        % Preliminary run - no control, 30 day build-up
        [Prelim, Prelim_ICs] = Get_ICs3(para0);
        
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
        [Classes] = ODEmodel3(para,Prelim_ICs); 

        % Post-Process for epidemic metrics
        [Classes, ~, Peak_hospital, ~, FinalHospital, ~, Days_lockdown, Days_Tier2, ~, nx, ix1, ix2, ~, ~] = PostProcessor3(Classes);
        Burden(idx,strat) = FinalHospital;
        Stringency(idx,strat) = ((para.LKRED^2)*Days_lockdown + (para.ICRED^2)*Days_Tier2);

        plotidx = 2*strat - 1;  % index for subfigure

        % plotting
        subplot(4,2,plotidx)

        yyaxis left
        ax1 = gca;
        ax1.YColor = 'k';
        ax1.FontSize = 16;
        ax1.FontSizeMode = 'manual';

        if idx == length(var_arr)
            plot(Classes.t, sum(Classes.IH,2), 'Color', 'r', 'Marker','none', 'LineWidth', 2.5, 'LineStyle','-')
            hold on
        else
            plot(Classes.t, sum(Classes.IH,2), 'Color', grayColor, 'Marker','none', 'LineWidth', 0.5, 'LineStyle','-')
            hold on
        end
    
        if plotidx > 6
            xlabel('Time (days)')
        end
        if plotidx == 1
            title('Active $I^H(t)$');
        end

        axis([0 maxtime 0 1250])

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

        axis([0 maxtime 0 1250])
        
        grid on

        plotidx = 2*strat;  % index for subfigure


        % plotting
        subplot(4,2,plotidx)

        if idx == length(var_arr)
            plot(Classes.t, sum(Classes.Hosp,2)./1000, 'Color', 'r', 'Marker','none', 'LineWidth', 2.5, 'LineStyle','-')
            hold on
        else
            plot(Classes.t, sum(Classes.Hosp,2)./1000, 'Color', grayColor, 'Marker','none', 'LineWidth', 0.5, 'LineStyle','-')
            hold on
        end
    
        if plotidx > 6
            xlabel('Time (days)')
        end
        if plotidx == 2
            title('Cumulative $I^H(t)$ (thousands)');
        end

        axis([0 maxtime 0 30])

        grid on

    end
end
toc

%save figure
saveas(gcf,strcat('./sens_images/','sim-sensitivity_',char(labs(var_to_test)),'.png'))


%% variability in cost

% strategy markers
markers = {'o','^','s','d'};
cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

% normalise cost contributions
Burden = Burden./max(Burden,[],2);
Stringency = Stringency./max(Stringency,[],2);

f = figure(2);
f.Position = [600 600 800 400];

for strat = strategies
    subplot(1,2,1)

    plot(var_arr(1:end-1),Burden(1:end-1,strat),"Color",'k','Marker',markers{strat},'MarkerSize',8, ...
         'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5)
    hold on
    axis([min(var_arr)*0.9 max(var_arr)*1.1 0 1.1])
    xlabel(strcat('$',vars(var_to_test),'$'))
    ylabel('Normalised cost contribution')
    title('Burden')

    if strat == length(strategies)
        xline(var_arr(end),'-','Model value','Color','r','Linewidth',2.5,'FontSize',16, ...
              'Interpreter','latex','LabelOrientation','horizontal','LabelVerticalAlignment','bottom')
    end


    subplot(1,2,2)

    plot(var_arr(1:end-1),Stringency(1:end-1,strat),"Color",'k','Marker',markers{strat},'MarkerSize',8, ...
         'MarkerFaceColor',cols(strat,:),'MarkerEdgeColor','k','LineWidth',1.5)
    hold on
    axis([min(var_arr)*0.9 max(var_arr)*1.1 0 1.1])
    xlabel(strcat('$',vars(var_to_test),'$'))
    title('Stringency')

    if strat == length(strategies)
        xline(var_arr(end),'-','Model value','Color','r','Linewidth',2.5,'FontSize',16, ...
              'Interpreter','latex','LabelOrientation','horizontal','LabelVerticalAlignment','bottom')
    end

end

%save figure
saveas(gcf,strcat('./sens_images/','cost-sensitivity_',char(labs(var_to_test)),'.png'))