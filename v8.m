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
delta = 0.06; % depriciation rate 
zeta = [0.95,0.93,0.89,0.83,0.73,0.57,0.38,0.21]; % survival rate 
zeta_inverse = [0.21,0.38,0.57,0.73,0.83,0.89,0.93,0.95];
% parameters
%  Capital grid
%  Be careful when setting maxkap across experiments! This value should not be binding! To see if it's binding, check the distribution of agents at k=maxkap in variables gkR and gkW (should be close to 0)
%  If you have problems with convergence, increasing nk might help.

maxkap = 2000000;                             % maximum value of capital grid  
minkap = 0;                             % minimum value of capital grid
nk=100;                                 % number of grid points
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
paramS.gamma_2 = -0.032;
paramS.gamma_1 = 0.32;
I_COST = 58000;
college_tuition = 12761;
elite_tuition = 33046;
alpha = 0.36; % output ela
elite_admit = [0.209,0.559,0.756];

paramS.mu_return = [0,0.25,0.47;0,0.20,0.35;0,0.28,0.56];
paramS.mu_high = [0,0,0];
paramS.mu_ec = [0.47,0.56,0.35];
paramS.mu_nc = [0.25,0.28,0.20];
paramS.mu_em = [0,0.25,0.47]; % high,oc,ec
paramS.mu_ib = [0,0.56,0.28]; %incorporated
paramS.mu_ub = [0,0.35,0.20];

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
k0=40000000;
l0=30000;

K0 =1;
r_guess=prod_em*alpha*K0^(alpha-1)-delta;
w_guess = prod_em*(1-alpha)*(K0^(alpha));
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
kapEduopt=ones(na,na,na,nk,nk,3);
neg=-1e10;  % very small number

fprintf('\nComputing general equilibrium...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital, labor and pension benefits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while q<nq && (abs(k1-k0)>tolk || abs(l1-l0)>tollab)
   
    q=q+1;
    
    fprintf('\nIteration %g out of %g \n',q,nq);
    
    % Prices
    r0  = prod_em*alpha*(k0^(alpha-1))*(l0^(1-alpha))-delta;
    w0 = prod_em*(1-alpha)*(k0^(alpha))*(l0^(-alpha));
    disp('  interest rate     wage');
    disp([r0, w0]);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Retired households    
    
    for i=J:-1:JR+1 %age
        for a_em=1:na
            for a_ub=1:na
                for a_ib =1:na
                    for e = 1:ne
                        for j=1:nk        % assets today
                            if i==J
                            % Initialize right-hand side of Bellman equation
                                kap0 =kap(j);
                                p = P * w0 * (a_em + exp(paramS.mu_em(e)+4 * paramS.gamma_1 + 22.7 * paramS.gamma_2));
                                cons1=(1+r0)*kap0+p;                                
                                income_n = prod_n*h_j(a_ib,e,2,paramS)*power(h_em(a_em,e,i,paramS),rho_ub)*kap0^nu_ub;
                                k_j_n = min((1+lambda_e)*kap0,power((delta+r0)/(income_n * nu_ub), 1/(nu_ub - 1)));                                
                                cons2 = (1-delta)*k_j_n+income_n*power(k_j_n,nu_ub)+(kap0-k_j_n)*(1+r0);
                                income_i = prod_i*h_j(a_ib,e,3,paramS)*power(h_em(a_em,e,i,paramS),rho_ib)*kap0^nu_ib;
                                k_j_i = min((1+lambda_e)*kap0,power((delta+r0)/(income_i * nu_ib), 1/(nu_ib - 1)));                                
                                cons3 = (1-delta)*k_j_i+income_i*power(k_j_i,nu_ib)+(kap0-k_j_i)*(1+r0);
                                if cons1<=0
                                    util1=neg;
                                else
                                    util1=(cons1^(1-sigma))/(1-sigma);
                                end    
                                if cons2<=0
                                    util2=neg;
                                else
                                    util2=(cons2^((1-sigma)))/(1-sigma);
                                end    
                                if cons3<=0
                                    util3=neg;
                                else
                                    util3=(cons3^((1-sigma)))/(1-sigma);
                                end    
                                vR(a_em,a_ub,a_ib,e,j,tR,1)=util1;
                                vR(a_em,a_ub,a_ib,e,j,tR,2)=util2;
                                vR(a_em,a_ub,a_ib,e,j,tR,3)=util3;                            
                            else
                                vmin=neg;
                                l=0;
                                % Loop over all k's in the capital grid to find the value,
                                % which gives max of the right-hand side of Bellman equation
                                while l<nk  	% assets tomorrow
                                    l=l+1;
                                    kap0=kap(j); % current asset holdings
                                    kap1=kap(l); % future asset holdings
                                    % Instantaneous utility
                                    p = P * w0 * (a_em + exp(paramS.mu_em(e)+4 * paramS.gamma_1 + 22.7 * paramS.gamma_2));
                                    cons=(1+r0)*kap0-kap1+p;
                                    
                                    if cons<=0
                                        util=neg;
                                    else
                                        util=(cons^((1-sigma)))/(1-sigma);
                                    end                                
    
                                    % Right-hand side of Bellman equation
                                    v0=util + beta*zeta(i-JR+1)*vR(a_em,a_ub,a_ib,e,l,i+1-JR,1); % (abi,abi,abi,edu,capital,time)
                                    
                                    % Store indirect utility and optimal saving 
                                    if v0>vmin
                                        vR(a_em,a_ub,a_ib,e,j,i-JR,1)=v0;
                                        kapRopt(a_em,a_ub,a_ib,e,j,i-JR,1)=l;
                                        vmin=v0;
                                    end
                                end
                                vmin=neg;
                                l=0;
                                while l<nk
                                    l=l+1;
                                    kap0=kap(j); % current asset holdings
                                    kap1=kap(l); % future asset holdings
                                    income_n = prod_n*h_j(a_ib,e,2,paramS)*power(h_em(a_em,e,i,paramS),rho_ub)*kap0^nu_ub;
                                    k_j_n = min((1+lambda_e)*kap0,power((delta+r0)/(income_n * nu_ub), 1/(nu_ub - 1)));                                
                                    cons = (1-delta)*k_j_n+income_n*power(k_j_n,nu_ub)+(kap0-k_j_n)*(1+r0)-kap1;
                                    if cons<=0
                                        util=neg;
                                    else
                                        util=(cons^((1-sigma)))/(1-sigma);
                                    end                                
                                    % Right-hand side of Bellman equation
                                    v0=util + zeta(i+1-JR)*beta*vR(a_em,a_ub,a_ib,e,l,i+1-JR,2);
                                    if v0>vmin
                                        vR(a_em,a_ub,a_ib,e,j,i-JR,2)=v0;
                                        kapRopt(a_em,a_ub,a_ib,e,j,i-JR,2)=l;
                                        vmin=v0;
                                    end
                                end
                                vmin=neg;
                                l=0;
                                while l<nk
                                    l=l+1;
                                    kap0=kap(j); % current asset holdings
                                    kap1=kap(l); % future asset holdings
    
                                    income_i = prod_i*h_j(a_ib,e,3,paramS)*power(h_em(a_em,e,i,paramS),rho_ib)*kap0^nu_ib;
                                    k_j_i = min((1+lambda_e)*kap0,power((delta+r0)/(income_i * nu_ib), 1/(nu_ib - 1)));                                
                                    cons = (1-delta)*k_j_i+income_i*power(k_j_i,nu_ib)+(kap0-k_j_i)*(1+r0)-kap1;
                                    if cons<=0
                                        util=neg;
                                    else
                                        util=(cons^((1-sigma)))/(1-sigma);
                                    end
                                    % Right-hand side of Bellman equation
                                    v0=util + zeta(i-JR)*beta*vR(a_em,a_ub,a_ib,e,l,i+1-JR,3);
                                    if v0>vmin
                                        vR(a_em,a_ub,a_ib,e,j,i-JR,3)=v0;
                                        kapRopt(a_em,a_ub,a_ib,e,j,i-JR,3)=l;
                                        vmin=v0;
                                    end
                                end
                            end
                            if vR(a_em,a_ub,a_ib,e,j,i-JR,3)>vR(a_em,a_ub,a_ib,e,j,i-JR,2) && vR(a_em,a_ub,a_ib,e,j,i-JR,3)>vR(a_em,a_ub,a_ib,e,j,i-JR,1)
                                careerR(a_em,a_ub,a_ib,e,j,i-JR) = 3;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i-JR) = vR(a_em,a_ub,a_ib,e,j,i-JR,3);
                            elseif vR(a_em,a_ub,a_ib,e,j,i-JR,2)>vR(a_em,a_ub,a_ib,e,j,i-JR,3) && vR(a_em,a_ub,a_ib,e,j,i-JR,2)>vR(a_em,a_ub,a_ib,e,j,i-JR,1)
                                careerR(a_em,a_ub,a_ib,e,j,i-JR) = 2;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i-JR) = vR(a_em,a_ub,a_ib,e,j,i-JR,2);
                            elseif vR(a_em,a_ub,a_ib,e,j,i-JR,1)>vR(a_em,a_ub,a_ib,e,j,i-JR,2) && vR(a_em,a_ub,a_ib,e,j,i-JR,1)>vR(a_em,a_ub,a_ib,e,j,i-JR,3)
                                careerR(a_em,a_ub,a_ib,e,j,i-JR) = 1;
                                vOptRetire(a_em,a_ub,a_ib,e,j,i-JR) = vR(a_em,a_ub,a_ib,e,j,i-JR,1);
                            end
                        end
                    end
                end
            end
        end       
    end

    % Working households
    for i=JR:-1:1               % age
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
                                cons=(1+r0)*kap0+w0*h_em(a_em,e,i,paramS)-kap1;                                
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                
                                % Right-hand side of Bellman equation
                                
                                % Need to be careful in terms of matrices, when transiting from last period of
                                % work to retirement:
                                
                                if i==JR       % retired next period
                                    v0=util + beta*vR(a_em,a_ub,a_ib,e,l,1,1);
                                else
                                    v0=util + beta*vW(a_em,a_ub,a_ib,e,l,i+1,1);
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
                                income_n = prod_n*h_j(a_ib,e,2,paramS)*power(h_em(a_em,e,i,paramS),rho_ub)*kap0^nu_ub;
                                k_j_n = min((1+lambda_e)*kap0,power((delta+r0)/(income_n * nu_ub), 1/(nu_ub - 1)));                                
                                cons = (1-delta)*k_j_n+income_n*power(k_j_n,nu_ub)+(kap0-k_j_n)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^(1-sigma))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                if i==JR       % retired next period
                                    v0=util + beta*vR(a_em,a_ub,a_ib,e,l,1,2);
                                else
                                    v0=util + beta*vW(a_em,a_ub,a_ib,e,l,i+1,2);
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

                                income_i = prod_i*h_j(a_ib,e,3,paramS)*power(h_em(a_em,e,i,paramS),rho_ib)*kap0^nu_ib;
                                k_j_i = min((1+lambda_e)*kap0,power((delta+r0)/(income_i * nu_ib), 1/(nu_ib - 1)));                                
                                cons = (1-delta)*k_j_i+income_i*power(k_j_i,nu_ib)+(kap0-k_j_i)*(1+r0)-kap1;
                                if cons<=0
                                    util=neg;
                                else
                                    util=(cons^((1-sigma)))/(1-sigma);
                                end
                                % Right-hand side of Bellman equation
                                if i==JR       % retired next period
                                    v0=util + beta*vR(a_em,a_ub,a_ib,e,l,1,3);
                                else
                                    v0=util + beta*vW(a_em,a_ub,a_ib,e,l,i+1,3);
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
                            else
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
                    for k=1:nk
                        v0=neg;
                        l=0;
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings     
                            % high school
                            cons = (1 + r0) * kap0-kap1;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^((1-sigma)))/(1-sigma);
                            end
                            v0= util + beta*vOptWork(a_em,a_ub,a_ib,1,l,1); % (abi,abi,abi,edu,capital,period)
                            if v0>vmin
                                % abi,abi,abi,edu,kap,choice
                                vEdu(a_em,a_ub,a_ib,j,k,1)=v0;
                                kapEduopt(a_em,a_ub,a_ib,j,k,1)=l;
                                vmin=v0;
                            end
                        % ordinary college
                        end

                        v0=neg;
                        l=0;
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings
                            k_fam=kap(k);
                            cons = (1 + r0) * (kap0 - college_tuition + financial_aid_nonelite(k_fam, a_em))-kap1;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^((1-sigma)))/(1-sigma);
                            end
                            % Right-hand side of Bellman equation
                            v0 = util + beta*vOptWork(a_em,a_ub,a_ib,2,l,1); % notice that the next period value function the capital is l!
                            if v0>vmin
                                vEdu(a_em,a_ub,a_ib,j,k,2)=v0;
                                kapEduopt(a_em,a_ub,a_ib,j,k,2)=l;
                                vmin=v0;
                            end
                        end
                        % elite college

                        v0=neg;
                        l=0;
                        while l<nk
                            l=l+1;
                            kap0=kap(j); % current asset holdings
                            kap1=kap(l); % future asset holdings
                            k_fam=kap(k);
                            cons = (1 + r0) * (kap0 -elite_tuition + financial_aid_elite(k_fam, a_em))-kap1;
                            if cons<=0
                                util=neg;
                            else
                                util=(cons^((1-sigma)))/(1-sigma);
                            end
                            % Right-hand side of Bellman equation
                            v0=util + beta*vOptWork(a_em,a_ub,a_ib,3,l,1);
                            if v0>vmin
                                vEdu(a_em,a_ub,a_ib,j,k,3)=v0;
                                kapEduopt(a_em,a_ub,a_ib,j,k,3)=l;
                                vmin=v0;
                            end
                        end
                        if vEdu(a_em,a_ub,a_ib,j,k,3)>vEdu(a_em,a_ub,a_ib,j,k,2) && vEdu(a_em,a_ub,a_ib,j,k,3)>vEdu(a_em,a_ub,a_ib,j,k,1)
                            careerEdu(a_em,a_ub,a_ib,j,k) = 3;
                            vOptEdu(a_em,a_ub,a_ib,j,k)=vEdu(a_em,a_ub,a_ib,j,k,3);
                        elseif vEdu(a_em,a_ub,a_ib,j,k,2)>vEdu(a_em,a_ub,a_ib,j,k,3) && vEdu(a_em,a_ub,a_ib,j,k,2)>vEdu(a_em,a_ub,a_ib,j,k,1)
                            careerEdu(a_em,a_ub,a_ib,j,k) = 2;
                            vOptEdu(a_em,a_ub,a_ib,j,k)=vEdu(a_em,a_ub,a_ib,j,k,2);
                        elseif vEdu(a_em,a_ub,a_ib,j,k,1)>vEdu(a_em,a_ub,a_ib,j,k,2) && vEdu(a_em,a_ub,a_ib,j,k,1)>vEdu(a_em,a_ub,a_ib,j,k,3)
                            careerEdu(a_em,a_ub,a_ib,j,k) = 1;
                            vOptEdu(a_em,a_ub,a_ib,j,k)=vEdu(a_em,a_ub,a_ib,j,k,1);
                        end
                    end
                end
            end
        end
    end
    %% Simulation %%
    fprintf("start simulation\n");
    [paramS.tranprob_abi_em,paramS.abi_loggrid_em,paramS.prob_abi_em] = markovappr(theta_em,abi_em_std,ra,na);
    [paramS.tranprob_abi_ub,paramS.abi_loggrid_ub,paramS.prob_abi_ub] = markovappr(theta_ub,abi_ub_std,ra,na);
    [paramS.tranprob_abi_ib,paramS.abi_loggrid_ib,paramS.prob_abi_ib] = markovappr(theta_ib,abi_ib_std,ra,na);    
    % simulation array
    nSim=5000;        
    BHistM = ones(nSim, J-2,2);
    cHistM = zeros(nSim, J);
    abiEMIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_em,paramS.prob_abi_em); % (index,age)
    abiUBIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_ub,paramS.prob_abi_ub);
    abiIBIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_ib,paramS.prob_abi_ib);
    k_idx_sim = round((nk-1).* rand(nSim,J)+1); % individual capital simulation ranges from [1,maxkap] size: nSim*J
    k_fam_idx_sim = round((nk-1).* rand(nSim,J)+1);
    KHistM = ones(nSim,J); % first period is k_idx_sim, we can calculate from schooling capital policy function.
    KHistM(:,1)=k_fam_idx_sim(:,1);
    eduSim =  ones(nSim,1);
    KEM =0; % held by employee
    BorrowSum=0;
    LHist = zeros(nSim,JR);
    % simulate first period: education decision
    fprintf("simulation first period\n");
    for t=1:J  
        for idxV = 1:nSim         
            a_em = abiEMIdx(idxV,t);
            a_ub = abiUBIdx(idxV,t);
            a_ib = abiIBIdx(idxV,t);
            kInd = k_idx_sim(idxV,t);
            k = k_fam_idx_sim(idxV,t);        
            edu_random_number = rand;         
            eduSim(idxV,t) = careerEdu(a_em,a_ub,a_ib,kInd,k);            
            % check if the student is admitted by elite college
            flag = 0;
            if (eduSim(idxV,t) == 3)
                sat_group=1;
	            if a_em >= 5
	                sat_group = 3;
                elseif a_em >= 2
	                sat_group = 2;
                end
	            if edu_random_number > elite_admit(sat_group)
		            eduSim(idxV,t) = careerEdu(a_em,a_ub,a_ib,kInd,k);
		            flag = 1; % admitted successfully.
                else
                    eduSim(idxV,t)=2;
                end            
            end
        end
    end   
    for t=1:J
        for idxV =1:nSim
            a_em = abiEMIdx(idxV,t);
            a_ub = abiUBIdx(idxV,t);
            a_ib = abiIBIdx(idxV,t);
            kInd = KHistM(idxV,t);            
            if t<=JR              
                W = careerW(a_em,a_ub,a_ib,e,kInd,t);
                KHistM(idxV,t+1)=kapWopt(a_em,a_ub,a_ib,e,kInd,t,W);
                if W==1
                    LHist(idxV,t)= h_em(a_em,e,t,paramS);                                        
                end
            else
                W = careerR(a_em,a_ub,a_ib,e,kInd,t-JR);
                KHistM(idxV,t+1)=kapRopt(a_em,a_ub,a_ib,e,KHistM(idxV,t),t-JR,W);            
            end            
        end
    end
    
   
    for t = 1 :JR    % age by age   
        for idxV =1:nSim
            e = eduSim(idxV,t);
            a_em = abiEMIdx(idxV,t);
            a_ub = abiUBIdx(idxV,t);
            a_ib = abiIBIdx(idxV,t);
            kInd = KHistM(idxV,t);
            kap0 = kap(kInd);
            W = careerW(a_em,a_ub,a_ib,e,kInd,t);    % why they do choose worker?         
            if W==1                
                KEM = KEM+kap0;
            elseif W==2                            
                income_n = prod_n*h_j(a_ib,e,2,paramS)*power(h_em(a_em,e,t,paramS),rho_ub)*kap0^nu_ub;
                k_j_n = min((1+kap0)*lambda_e,power((delta+r0)/(income_n * nu_ub), 1/(nu_ub - 1)));                
                BorrowSum = BorrowSum +k_j_n-kap0;  % at current period              
            elseif W==3                
                income_i = prod_i*h_j(a_ib,e,3,paramS)*power(h_em(a_em,e,t,paramS),rho_ib)*kap0^nu_ib;
                k_j_i = min((1+kap0)*lambda_e,power((delta+r0)/(income_i * nu_ib), 1/(nu_ib - 1)));
                BorrowSum = BorrowSum +k_j_i-kap0;                
            end
        end
    end
    for t=JR+1:J-1
        for idxV =1:nSim
            a_em = abiEMIdx(idxV,t);
            a_ub = abiUBIdx(idxV,t);
            a_ib = abiIBIdx(idxV,t);  
            e=eduSim(idxV,t);
            kInd = KHistM(idxV,t);
            kap0 = kap(kInd);
            W=careerR(a_em,a_ub,a_ib,e,kInd,t-JR);            
            if W==1
                KEM = KEM+kap0;
            elseif W==2            
                income_n = prod_n*h_j(a_ib,e,2,paramS)*power(h_em(a_em,e,t,paramS),rho_ub)*kap0^nu_ub;
                k_j_n = min((1+kap0)*lambda_e,power((delta+r0)/(income_n * nu_ub), 1/(nu_ub - 1)));                
                BorrowSum = BorrowSum +k_j_n-kap0;                    
            elseif W==3
                income_i = prod_i*h_j(a_ib,e,3,paramS)*power(h_em(a_em,e,t,paramS),rho_ib)*kap0^nu_ib;
                k_j_i = min((1+kap0)*lambda_e,power((delta+r0)/(income_i * nu_ib), 1/(nu_ib - 1)));
                BorrowSum = BorrowSum +k_j_i-kap0;                
            end
        end
    end
     % for id 
%% Update the guess on capital and labor
BSum = sum(BHistM,3);
k1 = max(0.01,KEM-BorrowSum);

k0=0.9*k0+0.1*k1;
l0=0.8*l0+0.2*mean(mean(LHist,1));

%% Display results
disp('  capital     labor');
disp([k0, l0]);
disp('deviation-capital deviation-labor       ');
disp([abs(k1-k0),  abs(l1-l0)]);
end 

%% Display equilibrium results
disp('      k0         l0       w         r    ');
disp([k0, l0, w0, r0]);

% Prices
r0  = prod_em*alpha*(k0^(alpha-1))*(l0^(1-alpha))-delta;
w0=prod_em*(1-alpha)*(k0^(alpha))*(l0^(-alpha));

% % Welfare
% welfare=0;
% 
% for i=1:tW
%     gkW_resh=reshape(gkW(:,i),N,nk);
%     sum_aux=gkW_resh.*vW(:,:,i);
%     welfare=welfare+sum(sum_aux(:));
% end
% 
% for i=1:tR
%     sum_aux=gkR.*vR;
%     welfare=welfare+sum(sum_aux(:));
% end
% 
% % Coefficient of variation    
% sum_aux=(kap-k0).*(kap-k0);
% sum_aux1=sum(gk');
% variance=sum(sum_aux1.*sum_aux);
% cv=sqrt(variance)/k0
% 
% %% Plots
% % Value function for a retired agent
% figure(1)
% age=5;
% plot1=plot(kap,vR(:,age));
% set(plot1,'LineWidth', 1.5);
% xlabel('capital','FontSize',14);
% ylabel('value function','FontSize',14);
% %title(['value function of a retired agent at age ', num2str(tW+age)])
% shg
% print('vf_retired','-dpng','-r300')
% 
% % Savings of a working agent
% figure(2)
% age=20;
% plot1=plot(kap,kap(kapWopt(1,:,age)),'k-',kap,kap(kapWopt(2,:,age)),'r--',kap,kap,'-x');
% set(plot1,'LineWidth', 1.0);
% xlabel('private asset holdings','FontSize',14);
% ylabel('saving','FontSize',14);
% legend('high shock', 'low shock','45 degree')
% %title(['saving of a working agent at age ', num2str(age)])
% shg
% print('saving_worker','-dpng','-r300')
% 
% % Savings profile
% figure(3)
% plot(1:J,kgen);
% xlabel('age','FontSize',14);
% ylabel('saving','FontSize',14);
% title('savings profile of workers')
% shg
% 
% % Labor supply of a working agent
% figure(5)
% age=20;
% plot(kap,labopt(1,:,age),'k-',kap,labopt(2,:,age),'r--');
% xlabel('capital','FontSize',14);
% ylabel('saving','FontSize',14);
% legend('high shock', 'low shock')
% title(['labor supply of a working agent at age ', num2str(age)])
% shg