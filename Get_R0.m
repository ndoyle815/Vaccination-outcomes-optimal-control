% function which uses the NextGeneration Matrix (NGM) approach to compute
% the basic reproduction number R_0 of the SEIR age model

function R = Get_R0(para)

% matrix elements ordered x = (E_{1-3}a, I_{1-3}a^A, I_{1-3}a^S, I_{1-3}a^PH, I_a^H)

T = zeros(13*para.n);
Sigma = zeros(13*para.n);

% transmission matrix T
T(1,3*para.n+1) = para.tau*para.beta(1,1);
T(1,3*para.n+2) = para.tau*para.beta(1,2);
T(1,3*para.n+3) = para.tau*para.beta(1,3);
T(1,4*para.n+1) = para.tau*para.beta(1,1);
T(1,4*para.n+2) = para.tau*para.beta(1,2);
T(1,4*para.n+3) = para.tau*para.beta(1,3);
T(1,5*para.n+1) = para.tau*para.beta(1,1);
T(1,5*para.n+2) = para.tau*para.beta(1,2);
T(1,5*para.n+3) = para.tau*para.beta(1,3);
T(1,6*para.n+1) = para.beta(1,1);
T(1,6*para.n+2) = para.beta(1,2);
T(1,6*para.n+3) = para.beta(1,3);
T(1,7*para.n+1) = para.beta(1,1);
T(1,7*para.n+2) = para.beta(1,2);
T(1,7*para.n+3) = para.beta(1,3);
T(1,8*para.n+1) = para.beta(1,1);
T(1,8*para.n+2) = para.beta(1,2);
T(1,8*para.n+3) = para.beta(1,3);
T(1,9*para.n+1) = para.beta(1,1);
T(1,9*para.n+2) = para.beta(1,2);
T(1,9*para.n+3) = para.beta(1,3);
T(1,10*para.n+1) = para.beta(1,1);
T(1,10*para.n+2) = para.beta(1,2);
T(1,10*para.n+3) = para.beta(1,3);
T(1,11*para.n+1) = para.beta(1,1);
T(1,11*para.n+2) = para.beta(1,2);
T(1,11*para.n+3) = para.beta(1,3);
T(1,12*para.n+1) = para.rho*para.beta(1,1);
T(1,12*para.n+2) = para.rho*para.beta(1,2);
T(1,12*para.n+3) = para.rho*para.beta(1,3);

T(2,3*para.n+1) = para.tau*para.beta(2,1);
T(2,3*para.n+2) = para.tau*para.beta(2,2);
T(2,3*para.n+3) = para.tau*para.beta(2,3);
T(2,4*para.n+1) = para.tau*para.beta(2,1);
T(2,4*para.n+2) = para.tau*para.beta(2,2);
T(2,4*para.n+3) = para.tau*para.beta(2,3);
T(2,5*para.n+1) = para.tau*para.beta(2,1);
T(2,5*para.n+2) = para.tau*para.beta(2,2);
T(2,5*para.n+3) = para.tau*para.beta(2,3);
T(2,6*para.n+1) = para.beta(2,1);
T(2,6*para.n+2) = para.beta(2,2);
T(2,6*para.n+3) = para.beta(2,3);
T(2,7*para.n+1) = para.beta(2,1);
T(2,7*para.n+2) = para.beta(2,2);
T(2,7*para.n+3) = para.beta(2,3);
T(2,8*para.n+1) = para.beta(2,1);
T(2,8*para.n+2) = para.beta(2,2);
T(2,8*para.n+3) = para.beta(2,3);
T(1,9*para.n+1) = para.beta(2,1);
T(1,9*para.n+2) = para.beta(2,2);
T(1,9*para.n+3) = para.beta(2,3);
T(1,10*para.n+1) = para.beta(2,1);
T(1,10*para.n+2) = para.beta(2,2);
T(1,10*para.n+3) = para.beta(2,3);
T(1,11*para.n+1) = para.beta(2,1);
T(1,11*para.n+2) = para.beta(2,2);
T(1,11*para.n+3) = para.beta(2,3);
T(1,12*para.n+1) = para.rho*para.beta(2,1);
T(1,12*para.n+2) = para.rho*para.beta(2,2);
T(1,12*para.n+3) = para.rho*para.beta(2,3);

T(3,3*para.n+1) = para.tau*para.beta(3,1);
T(3,3*para.n+2) = para.tau*para.beta(3,2);
T(3,3*para.n+3) = para.tau*para.beta(3,3);
T(3,4*para.n+1) = para.tau*para.beta(3,1);
T(3,4*para.n+2) = para.tau*para.beta(3,2);
T(3,4*para.n+3) = para.tau*para.beta(3,3);
T(3,5*para.n+1) = para.tau*para.beta(3,1);
T(3,5*para.n+2) = para.tau*para.beta(3,2);
T(3,5*para.n+3) = para.tau*para.beta(3,3);
T(3,6*para.n+1) = para.beta(3,1);
T(3,6*para.n+2) = para.beta(3,2);
T(3,6*para.n+3) = para.beta(3,3);
T(3,7*para.n+1) = para.beta(3,1);
T(3,7*para.n+2) = para.beta(3,2);
T(3,7*para.n+3) = para.beta(3,3);
T(3,8*para.n+1) = para.beta(3,1);
T(3,8*para.n+2) = para.beta(3,2);
T(3,8*para.n+3) = para.beta(3,3);
T(1,9*para.n+1) = para.beta(3,1);
T(1,9*para.n+2) = para.beta(3,2);
T(1,9*para.n+3) = para.beta(3,3);
T(1,10*para.n+1) = para.beta(3,1);
T(1,10*para.n+2) = para.beta(3,2);
T(1,10*para.n+3) = para.beta(3,3);
T(1,11*para.n+1) = para.beta(3,1);
T(1,11*para.n+2) = para.beta(3,2);
T(1,11*para.n+3) = para.beta(3,3);
T(1,12*para.n+1) = para.rho*para.beta(3,1);
T(1,12*para.n+2) = para.rho*para.beta(3,2);
T(1,12*para.n+3) = para.rho*para.beta(3,3);


% transition matrix Sigma
% E1,2,3
Sigma(1,1) = -para.n*para.epsilon;
Sigma(2,2) = -para.n*para.epsilon;
Sigma(3,3) = -para.n*para.epsilon;
Sigma(1*para.n+1,1) = para.n*para.epsilon;
Sigma(1*para.n+2,2) = para.n*para.epsilon;
Sigma(1*para.n+3,3) = para.n*para.epsilon;
Sigma(1*para.n+1,1*para.n+1) = -para.n*para.epsilon;
Sigma(1*para.n+2,1*para.n+2) = -para.n*para.epsilon;
Sigma(1*para.n+3,1*para.n+3) = -para.n*para.epsilon;
Sigma(2*para.n+1,1*para.n+1) = para.n*para.epsilon;
Sigma(2*para.n+2,1*para.n+2) = para.n*para.epsilon;
Sigma(2*para.n+3,1*para.n+3) = para.n*para.epsilon;
Sigma(2*para.n+1,2*para.n+1) = -para.n*para.epsilon;
Sigma(2*para.n+2,2*para.n+2) = -para.n*para.epsilon;
Sigma(2*para.n+3,2*para.n+3) = -para.n*para.epsilon;

% IA
Sigma(3*para.n+1,2*para.n+1) = (1 - para.da(1))*para.n*para.epsilon;
Sigma(3*para.n+2,2*para.n+2) = (1 - para.da(2))*para.n*para.epsilon;
Sigma(3*para.n+3,2*para.n+3) = (1 - para.da(3))*para.n*para.epsilon;
Sigma(3*para.n+1,3*para.n+1) = -para.n*para.gamma;
Sigma(3*para.n+2,3*para.n+2) = -para.n*para.gamma;
Sigma(3*para.n+3,3*para.n+3) = -para.n*para.gamma;
Sigma(4*para.n+1,3*para.n+1) = para.n*para.gamma;
Sigma(4*para.n+2,3*para.n+2) = para.n*para.gamma;
Sigma(4*para.n+3,3*para.n+3) = para.n*para.gamma;
Sigma(4*para.n+1,4*para.n+1) = -para.n*para.gamma;
Sigma(4*para.n+2,4*para.n+2) = -para.n*para.gamma;
Sigma(4*para.n+3,4*para.n+3) = -para.n*para.gamma;
Sigma(5*para.n+1,4*para.n+1) = para.n*para.gamma;
Sigma(5*para.n+2,4*para.n+2) = para.n*para.gamma;
Sigma(5*para.n+3,4*para.n+3) = para.n*para.gamma;
Sigma(5*para.n+1,5*para.n+1) = -para.n*para.gamma;
Sigma(5*para.n+2,5*para.n+2) = -para.n*para.gamma;
Sigma(5*para.n+3,5*para.n+3) = -para.n*para.gamma;


% IS
Sigma(6*para.n+1,2*para.n+1) = (1 - para.ha(1))*para.da(1)*para.n*para.epsilon;
Sigma(6*para.n+2,2*para.n+2) = (1 - para.ha(2))*para.da(2)*para.n*para.epsilon;
Sigma(6*para.n+3,2*para.n+3) = (1 - para.ha(3))*para.da(3)*para.n*para.epsilon;
Sigma(6*para.n+1,6*para.n+1) = -para.n*para.gamma;
Sigma(6*para.n+2,6*para.n+2) = -para.n*para.gamma;
Sigma(6*para.n+3,6*para.n+3) = -para.n*para.gamma;
Sigma(7*para.n+1,6*para.n+1) = para.n*para.gamma;
Sigma(7*para.n+2,6*para.n+2) = para.n*para.gamma;
Sigma(7*para.n+3,6*para.n+3) = para.n*para.gamma;
Sigma(7*para.n+1,7*para.n+1) = -para.n*para.gamma;
Sigma(7*para.n+2,7*para.n+2) = -para.n*para.gamma;
Sigma(7*para.n+3,7*para.n+3) = -para.n*para.gamma;
Sigma(8*para.n+1,7*para.n+1) = para.n*para.gamma;
Sigma(8*para.n+2,7*para.n+2) = para.n*para.gamma;
Sigma(8*para.n+3,7*para.n+3) = para.n*para.gamma;
Sigma(8*para.n+1,8*para.n+1) = -para.n*para.gamma;
Sigma(8*para.n+2,8*para.n+2) = -para.n*para.gamma;
Sigma(8*para.n+3,8*para.n+3) = -para.n*para.gamma;

% IPH
Sigma(9*para.n+1,2*para.n+1) = para.ha(1)*para.da(1)*para.n*para.epsilon;
Sigma(9*para.n+2,2*para.n+2) = para.ha(2)*para.da(2)*para.n*para.epsilon;
Sigma(9*para.n+3,2*para.n+3) = para.ha(3)*para.da(3)*para.n*para.epsilon;
Sigma(9*para.n+1,9*para.n+1) = -para.n*para.zeta;
Sigma(9*para.n+2,9*para.n+2) = -para.n*para.zeta;
Sigma(9*para.n+3,9*para.n+3) = -para.n*para.zeta;
Sigma(10*para.n+1,9*para.n+1) = para.n*para.zeta;
Sigma(10*para.n+2,9*para.n+2) = para.n*para.zeta;
Sigma(10*para.n+3,9*para.n+3) = para.n*para.zeta;
Sigma(10*para.n+1,10*para.n+1) = -para.n*para.zeta;
Sigma(10*para.n+2,10*para.n+2) = -para.n*para.zeta;
Sigma(10*para.n+3,10*para.n+3) = -para.n*para.zeta;
Sigma(11*para.n+1,10*para.n+1) = para.n*para.zeta;
Sigma(11*para.n+2,10*para.n+2) = para.n*para.zeta;
Sigma(11*para.n+3,10*para.n+3) = para.n*para.zeta;
Sigma(11*para.n+1,11*para.n+1) = -para.n*para.zeta;
Sigma(11*para.n+2,11*para.n+2) = -para.n*para.zeta;
Sigma(11*para.n+3,11*para.n+3) = -para.n*para.zeta;

% IH
Sigma(12*para.n+1,11*para.n+1) = para.n*para.zeta;
Sigma(12*para.n+2,11*para.n+2) = para.n*para.zeta;
Sigma(12*para.n+3,11*para.n+3) = para.n*para.zeta;
Sigma(12*para.n+1,12*para.n+1) = -para.delta;
Sigma(12*para.n+2,12*para.n+2) = -para.delta;
Sigma(12*para.n+3,12*para.n+3) = -para.delta;

% Next generation matrix
K = -T*inv(Sigma);
R = max(eig(K));