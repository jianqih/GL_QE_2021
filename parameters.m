% parameters
cS.num_id = 5000; % number of individual
cS.abi = 5; % grids
cS.shock = 5; % grids
cS.range_shock = 3;
cS.range_abi = 3;
cS.abi_em_std = exp(0.41);
cS.abi_ub_std = exp(0.32);
cS.abi_ib_std = exp(0.38);
cS.w_std = [0.66,0.73,0.59];
cS.w_h_std = 0.66;
cS.w_c_std = 0.73;
cS.w_e_std = 0.59;
cS.n_std = 0.59;
cS.i_std = 0.73;
% # ub:self-employment 0.0002
% # ib: entrepreneurship 0.0003
% # nc: ordinary college 0.020
% # ec: elite college 0.018
CONSUMP_N_STD = 0.0002; % non_incorp
CONSUMP_I_STD = 0.0003; % incorp
CONSUMP_C_STD = 0.020; % college
CONSUMP_E_STD = 0.018; % entrepreneur NC/EC
RAND_MAX = 1;
% K_CON_DIST = ones(5000,5000)*0.02;
CAREER_CHOICES = 3;
cS.theta_em = 0.47;
cS.theta_ib = 0.41;
cS.theta_ub = 0.38;
EDU = 3; %[no college, college, elite]
% ITER = 10;
cS.AGE = 17;
cS.RETIRE_AGE = 10; %5*9 = 45
EXP_BUS_TYPE_RETIRE = 3;
EXP_BUS_TYPE = 2;
cS.beta = 0.821; % discount rate
cS.delta = 0.266; % depriciation rate 
cS.zeta = [0.95,0.93,0.89,0.83,0.73,0.57,0.38,0.21]; % survival rate 
cS.sigma = 1.5;
cS.p = 0.4;
cS.nk = 9; % par.n
cS.kMin = -4; % left terminal of grid
cS.kMax = 4;
cS.kstep = 1;
cS.kGridV = linspace(cS.kMin, cS.kMax, cS.nk);  % Set k grid, it's a 1*nk row 
cS.kGridV = cS.kGridV';         % Make it to be a nk*1 vector
% K_FAM_STD = 0.5;
% COEF_ABI = ones([num_id,1])*0.3;
% COEF_K = 0.5; %
cS.lambda_e = 1.22; % constraint to entreur
NU_UB = 0.58; % lower than incorporated
NU_IB = 0.75;
RHO_UB = 0.03;
RHO_IB = 0.15;
PROD_N = 20.8;
PROD_I = 4.1;
PROD_EM = 2005;
ALPHA_2 = -0.032;
ALPHA_1 = 0.32;

I_COST = 58000;
COLLEGE_TUITION = 12761;
ELITE_TUITION = 33046;
alpha = 0.246; % output ela
CONSUMPTION_COLLEGE = 0;
CONSUMPTION_ELITE = 0;
% school choice
NO_COLLEGE = 1;
COLLEGE = 2;
ELITE = 3;
% job choice
WORKER = 1; 
NON_INCORP = 2;
INCORP = 3;
ELITE_ADMIT = [0.209,0.559,0.756];
cS.mu_high = [0,0,0];
cS.mu_ec = [0.47,0.56,0.35];
cS.mu_nc = [0.25,0.28,0.20];
cS.mu_em = [0,0.47,0.25];
cS.mu_ib = [0,0.56,0.28]; %incorporated
cS.mu_ub = [0,0.35,0.20];

[cS.tranprob_abi_em,cS.abi_grid_em,cS.prob_abi_em] = markovappr(cS.theta_em,cS.abi_em_std,cS.range_abi,cS.abi);
[cS.tranprob_abi_ub,cS.abi_grid_ub,cS.prob_abi_ub] = markovappr(cS.theta_ub,cS.abi_ub_std,cS.range_abi,cS.abi);
[cS.tranprob_abi_ib,cS.abi_grid_ib,cS.prob_abi_ib] = markovappr(cS.theta_ib,cS.abi_ib_std,cS.range_abi,cS.abi);

[~,cS.consump_grid_n,cS.prob_consump_n] = markovappr(0,CONSUMP_N_STD,cS.range_shock,cS.shock);
[~,cS.consump_grid_i,cS.prob_consump_i] = markovappr(0,CONSUMP_I_STD,cS.range_shock,cS.shock);
[~,cS.consump_grid_c,cS.prob_consump_c] = markovappr(0,CONSUMP_C_STD,cS.range_shock,cS.shock);
[~,cS.consump_grid_e,cS.prob_consump_e] = markovappr(0,CONSUMP_E_STD,cS.range_shock,cS.shock);

[~,cS.shock_grid_i,cS.prob_shock_i] = markovappr(0,cS.n_std,cS.range_shock,cS.shock);
[~,cS.shock_grid_n,cS.prob_shock_n] = markovappr(0,cS.i_std,cS.range_shock,cS.shock);

[~,cS.shock_grid_w1,cS.prob_shock_w1] = markovappr(0,cS.w_h_std,cS.range_shock,cS.shock);
[~,cS.shock_grid_w2,cS.prob_shock_w2] = markovappr(0,cS.w_c_std,cS.range_shock,cS.shock);
[~,cS.shock_grid_w3,cS.prob_shock_w3] = markovappr(0,cS.w_e_std,cS.range_shock,cS.shock);

cS.shock_grid_w = [cS.shock_grid_w1 cS.shock_grid_w2 cS.shock_grid_w3]';
cS.prob_shock_w = [cS.prob_shock_w1 cS.prob_shock_w2 cS.prob_shock_w3]';
cS.exp_em = [1:cS.AGE;0:cS.AGE-1;0:cS.AGE-1]';



h_em = [1.0000    1.6000    1.2840;
    2.0000    3.2000    2.5681;
    3.0000    4.8000    3.8521;
    4.0000    6.4000    5.1361;
    5.0000    8.0000    6.4201];

h_ib = [ 1.0000    1.7507    1.3231;
    2.0000    3.5013    2.6463;
    3.0000    5.2520    3.9694;
    4.0000    7.0027    5.2925;
    5.0000    8.7534    6.6156];


h_ub = [1.0000    1.4191    1.2214;
    2.0000    2.8381    2.4428;
    3.0000    4.2572    3.6642;
    4.0000    5.6763    4.8856;
    5.0000    7.0953    6.1070];


