% parameters
ABI = 5; % grids
SHOCK = 5; % grids
RANGE_SHOCK = 3;
RANGE_ABI = 3;
ABI_EM_STD = 0.38;
ABI_UB_STD = 0.38;
ABI_IB_STD = 0.32;
W_STD = [0.66,0.73,0.59];
W_H_STD = 0.66;
W_C_STD = 0.73;
W_E_STD = 0.59;
N_STD = 0.59;
I_STD = 0.73;
% # ub:self-employment 0.0002
% # ib: entrepreneurship 0.0003
% # nc: ordinary college 0.020
% # ec: elite college 0.018
CONSUMP_N_STD = 0.0002; % non_incorp
CONSUMP_I_STD = 0.0003; % incorp
CONSUMP_C_STD = 0.020; % college
CONSUMP_E_STD = 0.018; % entrepreneur NC/EC
RAND_MAX = 1;
K_CON_DIST = ones(5000,5000)*0.02;
CAREER_CHOICES = 3;
THETA_EM = 0.47;
THETA_IB = 0.41;
THETA_UB = 0.38;
EDU = 3; %[no college, college, elite]
ITER = 10;
AGE = 17;
RETIRE_AGE = 10; %5*9 = 45
EXP_BUS_TYPE_RETIRE = 3;
EXP_BUS_TYPE = 2;
beta = 0.821; % discount rate
delta = 0.266; % depriciation rate 
ZETA = [0.95,0.93,0.89,0.83,0.73,0.57,0.38,0.21]; % survival rate 
P = 0.4;
CAP = 20; % par.n
CAP_MIN = -4; % left terminal of grid
CAP_STEP = 1; % step
LAMBDA_E = 1.22; % constraint to entreur
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
ELITE_ADMIT = [0.209,0.559,0.756];

mu_high = [0,0,0];
mu_ec = [0.47,0.56,0.35];
mu_nc = [0.25,0.28,0.20];
mu_em = [0,0.47,0.25];
mu_ib = [0,0.56,0.28]; %incorporated
mu_ub = [0,0.35,0.20];

[tranprob_abi_em,abi_grid_em,prob_abi_em] = markovappr(THETA_EM,ABI_EM_STD,RANGE_ABI,ABI);
[tranprob_abi_ub,abi_grid_ub,prob_abi_ub] = markovappr(THETA_UB,ABI_UB_STD,RANGE_ABI,ABI);
[tranprob_abi_ib,abi_grid_ib,prob_abi_ib] = markovappr(THETA_IB,ABI_IB_STD,RANGE_ABI,ABI);

[~,consump_grid_n,prob_consump_n] = markovappr(0,CONSUMP_N_STD,RANGE_SHOCK,SHOCK);
[~,consump_grid_i,prob_consump_i] = markovappr(0,CONSUMP_I_STD,RANGE_SHOCK,SHOCK);
[~,consump_grid_c,prob_consump_c] = markovappr(0,CONSUMP_C_STD,RANGE_SHOCK,SHOCK);
[~,consump_grid_e,prob_consump_e] = markovappr(0,CONSUMP_E_STD,RANGE_SHOCK,SHOCK);

[~,shock_grid_n, prob_shock_i] = markovappr(0,N_STD,RANGE_SHOCK,SHOCK);
[~,shock_grid_i, prob_shock_n] = markovappr(0,I_STD,RANGE_SHOCK,SHOCK);

[~,shock_grid_w1,prob_shock_w1] = markovappr(0,W_STD(1),RANGE_SHOCK,SHOCK);
[~,shock_grid_w2,prob_shock_w2] = markovappr(0,W_STD(2),RANGE_SHOCK,SHOCK);
[~,shock_grid_w3,prob_shock_w3] = markovappr(0,W_STD(3),RANGE_SHOCK,SHOCK);

shock_grid_w = [shock_grid_w1;shock_grid_w2;shock_grid_w3];
prob_shock_w = [prob_shock_w1 prob_shock_w2 prob_shock_w3]';