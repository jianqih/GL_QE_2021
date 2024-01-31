% parameters
num_id = 5000; % number of individual
ABI = 5; % grids
SHOCK = 5; % grids
RANGE_SHOCK = 0.1; 
RANGE_ABI = 0.1;
ABI_EM_STD = exp(0.41);
ABI_UB_STD = exp(0.32);
ABI_IB_STD = exp(0.38);
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
BETA = 0.821; % discount rate
DELTA = 0.266; % depriciation rate 
ZETA = [0.95,0.93,0.89,0.83,0.73,0.57,0.38,0.21]; % survival rate 
SIGMA = 1.5;
P = 0.4;
CAP = 10; % par.n
CAP_MIN = 0; % left terminal of grid
CAP_MAX = 100000;
CAP_STEP = (CAP_MAX-CAP_MIN)/(CAP-1); % step
k=(CAP_MIN:CAP_STEP:CAP_MAX)';
K_FAM_STD = 0.5;
COEF_ABI = ones([num_id,1])*0.3;
COEF_K = 0.5; %
LAMBDA_E = 1.22; % constraint to entreur
NU_UB = 0.58; % lower than incorporated
NU_IB = 0.75;
RHO_UB = 0.03;
RHO_IB = 0.15;
PROD_N = 20.8;
PROD_I = 4.1;
PROD_EM = 2005;
ALPHA_2 = [-0.032,-0.032,-0.032];
ALPHA_1 = [0.32,0.32,0.32];
R = 0.05; %
W = 4000; % 

I_COST = [58000,8000];
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

mu_high = [0,0,0];
mu_ec = [0.47,0.56,0.35];
mu_nc = [0.25,0.28,0.20];
mu_em = [0,0.47,0.25];
mu_ib = [0,0.56,0.28]; %incorporated
mu_ub = [0,0.35,0.20];
