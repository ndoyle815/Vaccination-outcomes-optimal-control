% model parameters
gamma = 1/5.6;                  % infectious period
epsilon = 1/5.28;                % incubation period
omega = 1/750;                   % recovered period
zeta = gamma;                    % period between symptoms and hospital
delta = 1/8.78;                  % hospitalisation period
tau = 0.25;                      % relative infectiousness of asymptomatic
rho = 0.1;                       % relative infectiousness of hospitalised
da = [0.0275; 0.1350; 0.5374];   % probability of symptomatic infection
ha = [0.1268; 0.1173; 0.2430];   % probability of hospitalisation given symptomatic infection
N = [200000; 500000; 160000];    % population structure
n = size(N,1);                   % number of age classes
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity

% transmission matrix
beta = 0.945.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];

% initial exposed
E0 = 2e-4;

% gamma = 0.4641;
% zeta = gamma;
% beta = 2.45.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];
% E0 = 0.25e-4;
% N = 0.8.*N;

% gamma = 1/12;
% zeta = gamma;
% beta = 0.442.*[1.709, 0.458, 0.033; 0.497, 0.900, 0.073; 0.156, 0.374, 0.383];
% E0 = 10e-4;

% control parameters
strat = 1;                       % default strategy
init = 0;                        % default control state (no control)
tgap = 18;                       % decision-making gap (relaxation)
tdelay = 3;                      % natural delay in implementing decision
tdiff = 7;                       % tgap-tdiff = decision-making gap (reintroduction)
ICRED = 0.4;                     % reduction in transmission when in Intermediate Control
LKRED = 0.7;                     % reduction in transmission when in Lockdown

% default (dummy) switching thresholds
T01 = 10000;                     % No Control -> Intermediate Control
T10 = 10000;                     % Intermediate Control -> No Control
T12 = 10000;                     % Intermediate Control -> Lockdown
T21 = 10000;                     % Lockdown -> Intermediate Control

% Default time to run model for (preliminary run)
t0 = 0;
maxtime = 30;

% generate fixed vaccine parameters
% parameters varying sigmoid curve (ie. rollout speed)
tc = 100;                        % time to complete 50% vaccination
kappa = 0.05;                    % logistic shaping parameter
stagger = tc/2;                  % time between vaccination commencement of age groups
nu_a = [1;1;1];                  % vaccine transmission coefficient
efficacy = 0.9;                  % default vaccine efficacy
vstart = 2000;                   % default vaccine arrival date

save("./mats/Parameters.mat","beta","gamma","epsilon","omega","zeta","delta",...
    "tau","rho","da","ha","N","n","bedsper1000","Hmax","E0","strat","init",...
    "tgap","tdelay","tdiff","T01","T10","T12","T21","t0","maxtime","tc","kappa","stagger",...
    "nu_a","efficacy","vstart","ICRED","LKRED",'-mat')
