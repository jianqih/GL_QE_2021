% Instructor: Dean Corbae
% Fall 2016

% Macroeconomic consequences of eliminating social security in the U.S.
% This model is a simplified version of the model by Conesa and Krueger (1999). 

% This code is not necessarily written efficiently. It is written very
% explicitely, so that you can follow and understand it easier.

%% Clear the memory
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Demographics
J=17;                       % life-span
JR=10;                      % age of retirement
tR=J-JR+1;                  % length of retirement
tW=JR-1;                    % length of working life

% Preferences
sigma=1.5;                    % coefficient of relative risk aversion
beta = 0.821; % discount rate
delta = 0.266; % depriciation rate 
zeta = [0.95,0.93,0.89,0.83,0.73,0.57,0.38,0.21]; % survival rate 
zeta_inverse = [0.21,0.38,0.57,0.73,0.83,0.89,0.93,0.95];
% parameters
%  Capital grid
%  Be careful when setting maxkap across experiments! This value should not be binding! To see if it's binding, check the distribution of agents at k=maxkap in variables gkR and gkW (should be close to 0)
%  If you have problems with convergence, increasing nk might help.

maxkap = 10;                             % maximum value of capital grid  
minkap = 0;                             % minimum value of capital grid
nk=20;                                 % number of grid points
inckap=(maxkap-minkap)/(nk-1);       	% distance between points
aux=1:nk;
kap= minkap+inckap*(aux-1);  
na = 7;
ne = 3;
nb = 3;
abigrid = 1:na;
egrid = 1:ne;
[kg,ag,eg] = meshgrid(kap,abigrid,egrid);
ns = 5; % grids
rs = 3;
ra = 3;
abi_em_std = 0.38;
abi_ub_std = 0.38;
abi_ib_std = 0.32;
W_STD = [0.66,0.73,0.59];
W_H_STD = 0.66;
W_C_STD = 0.73;
W_E_STD = 0.59;
n_std = 0.59;
i_std = 0.73;
% # ub:self-employment 0.0002
% # ib: entrepreneurship 0.0003
% # nc: ordinary college 0.020
% # ec: elite college 0.018
CONSUMP_N_STD = 0.0002; % non_incorp
CONSUMP_I_STD = 0.0003; % incorp
CONSUMP_C_STD = 0.020; % college
CONSUMP_E_STD = 0.018; % entrepreneur NC/EC
RAND_MAX = 1;
CAREER_CHOICES = 3;
theta_em = 0.47;
theta_ib = 0.41;
theta_ub = 0.38;

P = 0.4;
lambda_e = 1.22; % constraint to entreur
nu_ub = 0.58; % lower than incorporated
nu_ib = 0.75;
rho_ub = 0.03;
rho_ib = 0.15;
prod_n = 20.8;
prod_i = 4.1;
prod_em = 2005;
gamma_2 = -0.032;
gamma_1 = 0.32;
I_COST = 58000;
college_tuition = 12761;
elite_tuition = 33046;
alpha = 0.246; % output ela
elite_admit = [0.209,0.559,0.756];

mu_high = [0,0,0];
mu_ec = [0.47,0.56,0.35];
mu_nc = [0.25,0.28,0.20];
mu_em = [0,0.47,0.25];
mu_ib = [0,0.56,0.28]; %incorporated
mu_ub = [0,0.35,0.20];

[tranprob_abi_em,abi_LogGrid_em,prob_abi_em] = markovappr(theta_em,abi_em_std,ra,na);
[tranprob_abi_ub,abi_LogGrid_ub,prob_abi_ub] = markovappr(theta_ub,abi_ub_std,ra,na);
[tranprob_abi_ib,abi_LogGrid_ib,prob_abi_ib] = markovappr(theta_ib,abi_ib_std,ra,na);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tolerance levels for capital, labor and pension benefits
tolk=1e-3;
tollab=1e-3;

nq=10;                                  % Max number of iterations
q=0;                                    % Counter for iterations

% Initial guesses for interest rate, wages and pension benefits
k0=3.4325;
l0=0.3263;

k1=k0+10;
l1=l0+10;

% Initializations for backward induction
vR=zeros(na,na,na,ne,nk,tR,3);                        % value function of retired agents
vOptRetire=zeros(na,na,na,ne,nk,tR);
kapRopt=ones(na,na,na,ne,nk,tR,3);                    % optimal savings of retired agents
% (store INDEX of k' in capital grid, not k' itself!)
careerR=ones(na,na,na,ne,nk,tW,3);
careerW=ones(na,na,na,ne,nk,tW);
careerEdu=ones(na,na,na,nk,nk,3);
vW=zeros(na,na,na,ne,nk,tW,3);                      % value function of workers
vOptWork=zeros(na,na,na,ne,nk,tW);                      % value function of workers

kapWopt=ones(na,na,na,ne,nk,tW,3);                  % optimal savings of workers

% (store INDEX of k' in capital grid, not k' itself!)
vEdu=zeros(na,na,na,nk,nk,3);
vOptEdu=zeros(na,na,na,nk,nk);
kapEduOpt=ones(na,na,na,nk,nk,3);
neg=-1e10;  % very small number

fprintf('\nComputing general equilibrium...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital, labor and pension benefits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while q<nq && (abs(k1-k0)>tolk || abs(l1-l0)>tollab)
   
    q=q+1;
    
    fprintf('\nIteration %g out of %g \n',q,nq);
    
    % Prices
    r0  = alpha*(k0^(alpha-1))*(l0^(1-alpha))-delta;
    w0 = (1-alpha)*(k0^(alpha))*(l0^(-alpha));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Retired households
   
    % Last period utility    
    cons=(1+r0)*kg(:) + P*w0*(ag(:)+exp(mu_em(eg(:))')) * exp(4 * gamma_1 + 22.7 * gamma_2);                    % last period consumption (vector!)
    util=(cons.^(1-sigma))/(1-sigma);       % last period utility (vector!)
    util_exp=repmat(util,[na*na,1]);
    vR(:,:,:,:,:,tR)=reshape(util_exp,na,na,na,ne,nk);                       % last period indirect utility (vector!)
    
    for i=tR-1:-1:1       % age
        for a_em=1:na
            for a_ub=1:na
                for a_ib =1:na
                    for e = 1:ne
                        for j=1:nk        % assets today
                           
                            % Initialize right-hand side of Bellman equation
                            vmin=neg;
                            l=0;
                            
                            % Loop over all k's in the capital grid to find the value,
                            % which gives max of the right-hand side of Bellman equation
                            
                            while l<nk  	% assets tomorrow
                                l=l+1;
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings
                                
                                % Instantaneous utility
                                p = P * w0 * (a_em + exp(mu_em(e)+4 * gamma_1 + 22.7 * gamma_2));
                                cons=(1+r0)*kap0-kap1+p;
                                
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                
                                % Right-hand side of Bellman equation
                                v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,i+1);
                                
                                % Store indirect utility and optimal saving 
                                if v0>vmin
                                    vR(a_em,a_ub,a_ib,e,j,i,1)=v0;
                                    kapRopt(a_em,a_ub,a_ib,e,j,i)=l;
                                    vmin=v0;
                                end
                            end
                            vmin=neg;
                            l=0;
                            while l<nk
                                l=l+1;
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings

                                i_scaler_non_incorp = prod_n*(a_ub+exp(mu_ub(e)))*power(a_em + exp(mu_em(e)+exp_em(i,e)*gamma_1 + exp_em(i,e)*gamma_2*exp_em(i,e)),rho_ub)*kap0^nu_ub;
                                invest_non_incorp = min((1+lambda_e)*kap0,power((delta+r0)/(i_scaler_non_incorp * nu_ub), 1/(nu_ub - 1)));                                
                                cons = (1-delta)*invest_non_incorp+i_scaler_non_incorp*power(invest_non_incorp,nu_ub)+(kap0-invest_non_incorp)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,i+1);
                                if v0>vmin
                                    vR(a_em,a_ub,a_ib,e,j,i,2)=v0;
                                    kapRopt2(a_em,a_ub,a_ib,e,j,i)=l;
                                    vmin=v0;
                                end
                            end
                            vmin=neg;
                            l=0;
                            while l<nk
                                l=l+1;
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings

                                i_scaler_incorp = prod_i*(a_ib+exp(mu_ib(e)))*power(a_em + exp(mu_em(e)+exp_em(i,e)*gamma_1 + exp_em(i,e)*gamma_2*exp_em(i,e)),rho_ib)*kap0^nu_ib;
                                invest_incorp = min((1+lambda_e)*kap0,power((delta+r0)/(i_scaler_incorp * nu_ib), 1/(nu_ib - 1)));                                
                                cons = (1-delta)*invest_incorp+i_scaler_incorp*power(invest_incorp,nu_ib)+(kap0-invest_incorp)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,i+1);
                                if v0>vmin
                                    vR(a_em,a_ub,a_ib,e,j,i,3)=v0;
                                    kapRopt3(a_em,a_ub,a_ib,e,j,i)=l;
                                    vmin=v0;
                                end
                            end
                            if vR(a_em,a_ub,a_ib,e,j,i,3)>vR(a_em,a_ub,a_ib,e,j,i,2) && vR(a_em,a_ub,a_ib,e,j,i,3)>vR(a_em,a_ub,a_ib,e,j,i,1)
                                careerR(a_em,a_ub,a_ib,e,j,i) = 3;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i) = vR(a_em,a_ub,a_ib,e,j,i,3);
                            elseif vR(a_em,a_ub,a_ib,e,j,i,2)>vR(a_em,a_ub,a_ib,e,j,i,3) && vR(a_em,a_ub,a_ib,e,j,i,2)>vR(a_em,a_ub,a_ib,e,j,i,1)
                                careerR(a_em,a_ub,a_ib,e,j,i) = 2;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i) = vR(a_em,a_ub,a_ib,e,j,i,2);
                            elseif vR(a_em,a_ub,a_ib,e,j,i,1)>vR(a_em,a_ub,a_ib,e,j,i,2) && vR(a_em,a_ub,a_ib,e,j,i,1)>vR(a_em,a_ub,a_ib,e,j,i,3)
                                careerR(a_em,a_ub,a_ib,e,j,i) = 1;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i) = vR(a_em,a_ub,a_ib,e,j,i,1);
                            end
                        end
                    end
                end
            end
        end       
    end

    % Working households
    for i=tW:-1:1               % age
        for a_em=1:na
            for a_ub=1:na
                for a_ib =1:na
                    for e = 1:ne
                        for j=1:nk        % assets today                
                            % Initialize right-hand side of Bellman equation
                            vmin=neg;
                            l=0;
                            
                            % Loop over all k's in the capital grid to find the value,
                            % which gives max of the right-hand side of Bellman equation
                            
                            while l<nk  	% assets tomorrow
                                l=l+1;
                                
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings                                                               
                                
                                % Instantaneous utility
                                cons=(1+r0)*kap0+w0*(a_em + exp(mu_em(e)+exp_em(i,e)*gamma_1+gamma_2*exp_em(i,e)*exp_em(i,e)))-kap1;
                                
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^(1-sigma))/(1-sigma);
                                end
                                
                                % Right-hand side of Bellman equation
                                
                                % Need to be careful in terms of matrices, when transiting from last period of
                                % work to retirement:
                                
                                if i==tW       % retired next period
                                    v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,1);
                                else
                                    v0=util + beta*vOptWork(a_em,a_ub,a_ib,e,j,i+1);
                                end
                                
                                % Store indirect utility, optimal saving and labor
                                if v0>vmin
                                    vW(a_em,a_ub,a_ib,e,j,i,1)=v0;
                                    kapWopt(a_em,a_ub,a_ib,e,j,i,1)=l;                            
                                    vmin=v0;
                                end
                            end
                            vmin=neg;
                            l=0;
                            while l<nk
                                l=l+1;
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings
                                i_scaler_non_incorp = prod_n*(a_ub+exp(mu_ub(e)))*power(a_em + exp(mu_em(e)+exp_em(i,e)*gamma_1 + exp_em(i,e)*gamma_2*exp_em(i,e)),rho_ub)*kap0^nu_ub;
                                invest_non_incorp = min((1+lambda_e)*kap0,power((delta+r0)/(i_scaler_non_incorp * nu_ub), 1/(nu_ub - 1)));                                
                                cons = (1-delta)*invest_non_incorp+i_scaler_non_incorp*power(invest_non_incorp,nu_ub)+(kap0-invest_non_incorp)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                if i==tW       % retired next period
                                    v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,1);
                                else
                                    v0=util + beta*vOptWork(a_em,a_ub,a_ib,e,j,i+1);
                                end
                                if v0>vmin
                                    vW(a_em,a_ub,a_ib,e,j,i,2)=v0;
                                    kapWopt(a_em,a_ub,a_ib,e,j,i,2)=l;
                                    vmin=v0;
                                end
                            end
                            vmin=neg;
                            l=0;
                            while l<nk
                                l=l+1;
                                kap0=kap(j); % current asset holdings
                                kap1=kap(l); % future asset holdings

                                i_scaler_incorp = prod_i*(a_ib+exp(mu_ib(e)))*power(a_em + exp(mu_em(e)+exp_em(i,e)*gamma_1 + exp_em(i,e)*gamma_2*exp_em(i,e)),rho_ib)*kap0^nu_ib;
                                invest_incorp = min((1+lambda_e)*kap0,power((delta+r0)/(i_scaler_incorp * nu_ib), 1/(nu_ib - 1)));                                
                                cons = (1-delta)*invest_incorp+i_scaler_incorp*power(invest_incorp,nu_ib)+(kap0-invest_incorp)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                if i==tW       % retired next period
                                    v0=util + beta*vOptRetire(a_em,a_ub,a_ib,e,j,1);
                                else
                                    v0=util + beta*vOptWork(a_em,a_ub,a_ib,e,j,i+1);
                                end
                                if v0>vmin
                                    vW(a_em,a_ub,a_ib,e,j,i,3)=v0;
                                    kapWopt(a_em,a_ub,a_ib,e,j,i,3)=l;
                                    vmin=v0;
                                end
                            end
                            if vW(a_em,a_ub,a_ib,e,j,i,3)>vW(a_em,a_ub,a_ib,e,j,i,2) && vW(a_em,a_ub,a_ib,e,j,i,3)>vW(a_em,a_ub,a_ib,e,j,i,1)
                                careerW(a_em,a_ub,a_ib,e,j,i) = 3;
                                vOptWork(a_em,a_ub,a_ib,e,j,i)=vW(a_em,a_ub,a_ib,e,j,i,3);
                            elseif vW(a_em,a_ub,a_ib,e,j,i,2)>vW(a_em,a_ub,a_ib,e,j,i,3) && vW(a_em,a_ub,a_ib,e,j,i,2)>vW(a_em,a_ub,a_ib,e,j,i,1)
                                careerW(a_em,a_ub,a_ib,e,j,i) = 2;
                                vOptWork(a_em,a_ub,a_ib,e,j,i)=vW(a_em,a_ub,a_ib,e,j,i,2);
                            elseif vW(a_em,a_ub,a_ib,e,j,i,1)>vW(a_em,a_ub,a_ib,e,j,i,2) && vW(a_em,a_ub,a_ib,e,j,i,1)>vW(a_em,a_ub,a_ib,e,j,i,3)
                                careerW(a_em,a_ub,a_ib,e,j,i) = 1;
                                vOptWork(a_em,a_ub,a_ib,e,j,i)=vW(a_em,a_ub,a_ib,e,j,i,1);
                            end
                        end
                    end
                end
            end
        end
    end
    
    fprintf("start schooling decision\n");
    for a_em=1:na
        for a_ub=1:na
            for a_ib =1:na
                for j = 1:nk
                    for k_fam=1:nk
                        v0=neg;
                        l=0;
                        
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings
                            cons = (1 + r0) * kap0;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^((1-sigma)))/(1-sigma);
                            end
                            v0 = util+beta*vOptWork(a_em,a_ub,a_ib,1,j,1);
                            if v0>vmin
                                vEdu(a_em,a_ub,a_ib,j,k_fam,1)=v0;
                                kapEduOpt(a_em,a_ub,a_ib,j,k_fam,1)=l;
                                vmin=v0;
                            end
                        end
                        v0=neg;
                        l=0;
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings
                            cons = (1 + r0) * (kap0 - 0.0001*college_tuition + financial_aid_nonelite(k_fam, a_em))-kap1;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^(1-sigma))/(1-sigma);
                            end
                            % Right-hand side of Bellman equation
                            v0=util + beta*vOptWork(a_em,a_ub,a_ib,2,j,1);
                            if v0>vmin
                                vEdu(a_em,a_ub,a_ib,j,k_fam,2)=v0;
                                kapEduOpt(a_em,a_ub,a_ib,j,k_fam,2)=l;
                                vmin=v0;
                            end
                        end
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings
                            cons = (1 + r0) * (kap0 - 0.0001*elite_tuition + financial_aid_nonelite(k_fam, a_em))-kap1;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^(1-sigma))/(1-sigma);
                            end
                            % Right-hand side of Bellman equation
                            v0=util + beta*vW(a_em,a_ub,a_ib,3,j,1);
                            if v0>vmin
                                vEdu(a_em,a_ub,a_ib,j,k_fam,3)=v0;
                                kapEduOpt(a_em,a_ub,a_ib,j,k_fam,3)=l;
                                vmin=v0;
                            end
                        end
                        if vEdu(a_em,a_ub,a_ib,j,k_fam,3)>vEdu(a_em,a_ub,a_ib,j,k_fam,2) && vEdu(a_em,a_ub,a_ib,j,k_fam,3)>vEdu(a_em,a_ub,a_ib,j,k_fam,1)
                            careerEdu(a_em,a_ub,a_ib,j,k_fam) = 3;
                            vOptEdu(a_em,a_ub,a_ib,j,k_fam)=vEdu(a_em,a_ub,a_ib,j,k_fam,3);
                        elseif vEdu(a_em,a_ub,a_ib,j,k_fam,2)>vEdu(a_em,a_ub,a_ib,j,k_fam,3) && vEdu(a_em,a_ub,a_ib,j,k_fam,2)>vEdu(a_em,a_ub,a_ib,j,k_fam,1)
                            careerEdu(a_em,a_ub,a_ib,j,k_fam) = 2;
                            vOptEdu(a_em,a_ub,a_ib,j,k_fam)=vEdu(a_em,a_ub,a_ib,j,k_fam,2);
                        elseif vEdu(a_em,a_ub,a_ib,j,k_fam,1)>vEdu(a_em,a_ub,a_ib,j,k_fam,2) && vEdu(a_em,a_ub,a_ib,j,k_fam,1)>vEdu(a_em,a_ub,a_ib,j,k_fam,3)
                            careerEdu(a_em,a_ub,a_ib,j,k_fam) = 1;
                            vOptEdu(a_em,a_ub,a_ib,j,k_fam)=vEdu(a_em,a_ub,a_ib,j,k_fam,1);
                        end
                    end
                end
            end
        end
    end
end
