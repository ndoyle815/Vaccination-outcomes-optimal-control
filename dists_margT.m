% Script to evaluate control strategies with vaccination, taking the
% arrival date distributions into account and producing Figure 5+6 from the
% report
clear all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination start times
vstart_times = [180:10:1080];

% Define time to run model for
t_init = 30;     % preliminary run
maxtime = max(vstart_times);  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 100 150 600; 50 100 150 200; 100 300 500 600; 250 350 425 500];
strategies = [1:length(thresholds)];

% add control thresholds defined by strategy
para = para0;
para.maxtime = maxtime;
para.Hmax = 1500;        % modify hospital capacity

% define vaccine efficacy
para.efficacy = 0.1;

% define functional weights
weights = [0.0:0.01:1.0];
w2 = 2;

% stores cost function outputs
ns = length(strategies);
nw = length(weights);
nv = length(vstart_times);
fs = zeros(nw,nv,ns);

tic
for strat = strategies
    strat
    % set switching thresholds
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % run preliminary simulation to get ICs
    [Prelim, Prelim_ICs] = Get_ICs(para0);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.T12
        para.init = 1;
    else
        para.init = 2;
    end

    for v = 1:nv
        para.vstart = vstart_times(v);

        % Run main simulation
        [Classes, burden, stringency, peak_hospital] = ODEmodel(para, Prelim_ICs);
    
        for w = 1:nw
            % evaluate cost function
            fs(w,v,strat) = CostFunction([weights(w), w2], para, burden, stringency, peak_hospital, 0);
        end
    end
end
toc

%% Plotting
set(0,'defaultaxesfontsize',20)

f = figure(1);
f.Position = [100 400 1200 400];

cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

subplot(1,3,[1 2])
ax = gca;
ax.Position(1) = ax.Position(1) - 0.05;
ax.Position(2) = ax.Position(2) + 0.02;

for strat = strategies
    fstrat = fs(:,:,strat);
    h1 = surf(fstrat,'FaceColor',cols(strat,:),'LineStyle','none');
    hold on
    
    %%Extract X,Y and Z data from surface plot
    x=h1.XData;
    y=h1.YData;
    z=h1.ZData;
    %%Create vectors out of surface's XData and YData
    x=x(1,:);
    y=y(:,1);
    %%Divide the lengths by the number of lines needed
    xnumlines = ceil(length(vstart_times)/5); % 10 lines
    ynumlines = ceil(length(weights)/5); % 10 partitions
    xspacing = round(length(x)/xnumlines);
    yspacing = round(length(y)/ynumlines);
    %%Plot the mesh lines 
    % Plotting lines in the X-Z plane
    for i = 1:yspacing:length(y)
        Y1 = y(i)*ones(size(x)); % a constant vector
        Z1 = z(i,:);
        plot3(x,Y1,Z1,'-k','LineWidth',0.5);
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(x)
        X2 = x(i)*ones(size(y)); % a constant vector
        Z2 = z(:,i);
        plot3(X2,y,Z2,'-k','LineWidth',0.5);
    end
end
view([315 45])
yticks([1:20:nw])
xticks([1:18:nv])
xtickangle(0)
ytickangle(0)
yticklabels(weights(1:20:end-5))
xticklabels(vstart_times(1:18:end))
ylim([1 nw])
xlim([1 nv])
zticks([0.4 0.6 0.8 1])
zlim([0.3 1])
label_y = ylabel('Weight $w$','HorizontalAlignment','right','Rotation',0,'VerticalAlignment','baseline');
label_x = xlabel('$T$ (days)','HorizontalAlignment','left','Rotation',0,'VerticalAlignment','baseline');
f.CurrentAxes.ZDir = 'Reverse';
zlabel('Cost')

set(0,'defaultaxesfontsize',21.6)

subplot(1,3,3)
ax = gca;
%ax.Position(1) = ax.Position(1) + 0.03;
ax.Position(2) = ax.Position(2) + 0.12;
ax.Position(4) = ax.Position(4) - 0.12;

for strat = strategies
    fstrat = fs(:,:,strat);
    h2(strat) = surf(fstrat,'FaceColor','none','FaceColor',cols(strat,:),'LineStyle','none');
    hold on

    %%Extract X,Y and Z data from surface plot
    h2strat = h2(strat);
    x=h2strat.XData;
    y=h2strat.YData;
    z=h2strat.ZData;
    %%Create vectors out of surface's XData and YData
    x=x(1,:);
    y=y(:,1);
    %%Divide the lengths by the number of lines needed
    xnumlines = ceil(length(vstart_times)/5); % 10 lines
    ynumlines = ceil(length(weights)/5); % 10 partitions
    xspacing = round(length(x)/xnumlines);
    yspacing = round(length(y)/ynumlines);
    %%Plot the mesh lines 
    % Plotting lines in the X-Z plane
    for i = 1:yspacing:length(y)
        Y1 = y(i)*ones(size(x)); % a constant vector
        Z1 = z(i,:);
        plot3(x,Y1,Z1,'-k','LineWidth',0.5);
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(x)
        X2 = x(i)*ones(size(y)); % a constant vector
        Z2 = z(:,i);
        plot3(X2,y,Z2,'-k','LineWidth',0.5);
    end
end
view([0 270])
f.CurrentAxes.YDir = 'Reverse';
yticks([1:20:nw])
xticks([1:18:nv])
xtickangle(0)
yticklabels(weights(1:20:end))
xticklabels(vstart_times(1:18:end))
ylim([1 nw])
xlim([1 nv])
ylabel('Weight $w$','FontSize',22)
label_x = xlabel('$T$ (days)','FontSize',22);
label_x.Position(2) = label_x.Position(2) - 0.6;
if para.efficacy == 0.3
    legend(h2, 'S1 (Cautious easing)', 'S2 (Suppression)', 'S3 (Slow control)', 'S4 (Rapid control)','Location','eastoutside','Interpreter','latex')
end

saveas(f,strcat('./vacc_images/optstrat_eta',num2str(para.efficacy),'.png'));


