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
J=66;                       % life-span
JR=46;                      % age of retirement
tR=J-JR+1;                  % length of retirement
tW=JR-1;                    % length of working life
n=0.011;                    % Population growth

% Preferences
beta=0.97;                  % discount factor
sigma=2;                    % coefficient of relative risk aversion
gamma=0.42;                 % weight on consumption

% Production
alpha=0.36;                 % production elasticity of capital 
delta=0.06;                 % rate of depreciation

p = 0.4;            % pension
% Social Security tax rate
theta=0.11;

% Idiosyncratic productivities
N=2;                        % number of shock realizations at each age

% Age-efficiency profile
eff=load('ef.txt');

% Distribution of newborns over shocks
z_init(1)=0.2037;
z_init(2)=0.7963;

% Idiosyncratic productivities
z=[3.0, 0.5];

% Transition matrix
tran(1,1)=0.9261;
tran(1,2)=1.0-0.9261;
tran(2,2)=0.9811;
tran(2,1)=1.0-0.9811;

% Measure of each generation
mass=ones(J,1);
for j=2:J
    mass(j)=mass(j-1)/(1+n);
end

% Normalized measure of each generation (sum up to 1)
mass=mass/sum(mass);

%  Capital grid
%  Be careful when setting maxkap across experiments! This value should not be binding! To see if it's binding, check the distribution of agents at k=maxkap in variables gkR and gkW (should be close to 0)
%  If you have problems with convergence, increasing nk might help.

maxkap =20;                             % maximum value of capital grid  
minkap = 0;                             % minimum value of capital grid
nk=200;                                 % number of grid points
inckap=(maxkap-minkap)/(nk-1);       	% distance between points
aux=1:nk;
kap= minkap+inckap*(aux-1);  
na = 5;
ne = 3;
nb = 3;
abigrid = 1:na;
egrid = 1:ne;
[kg,ag1,ag2,ag3,eg] = meshgrid(kap,abigrid,abigrid,abigrid,egrid);


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
vR=zeros(nk,tR);                        % value function of retired agents
kapRopt=ones(nk,tR);                    % optimal savings of retired agents
% (store INDEX of k' in capital grid, not k' itself!)

vW=zeros(N,nk,tW);                      % value function of workers
kapWopt=ones(N,nk,tW);                  % optimal savings of workers
% (store INDEX of k' in capital grid, not k' itself!)

labopt=ones(N,nk,tW);                   % optimal labor supply 

neg=-1e10;                              % very small number

fprintf('\nComputing general equilibrium...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over capital, labor and pension benefits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while q<nq && (abs(k1-k0)>tolk || abs(l1-l0)>tollab)
   
    q=q+1;
    
    fprintf('\nIteration %g out of %g \n',q,nq);
    
    % Prices
    r0  = alpha*(k0^(alpha-1))*(l0^(1-alpha))-delta;
    w0=(1-alpha)*(k0^(alpha))*(l0^(-alpha));
    
    % Pension benefit
    % ss0=theta*w0*l0/sum(mass(JR:end));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% BACKWARD INDUCTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Retired households
   
    
    
    for i=tR:-1:1       % age
        for e = 1:ne
            for a_em = 1:na
                for a_ub = 1:na
                    for a_ib = 1:na % abi ib                       
                        for j = 1:nk        % assets today
                            for shock_n = 1:N
                                for shock_i = 1:N
                                    h_em_ = h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e)));
                                    i_scaler_non_incorp = prod(e) * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_,rho);
                                    i_scaler_incorp = prod(e) * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, rho);
                                    invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp .* NU_UB), 1/(NU_UB - 1));
                                    k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) + (kap(k) - invest_non_incorp) * (1 + R);
                                    if i == tR
                                    % Last period utility
                                        cons1=(1+r0)*kg(:) + P*w0*h_em(ag(:),eg(:)) * exp((4 * ALPHA_1 + 22.7 * ALPHA_2));   % last period consumption (vector!)
                                        cons2=(1+r0)*kg(:) + P*w0*h_em(ag(:),eg(:)) * exp((4 * ALPHA_1 + 22.7 * ALPHA_2));   % last period consumption (vector!)
                                        cons3=(1+r0)*kg(:) + P*w0*h_em(ag(:),eg(:)) * exp((4 * ALPHA_1 + 22.7 * ALPHA_2));   % last period consumption (vector!)
                                        util=(cons.^(1-sigma))/(1-sigma);       % last period utility (vector!)
                                        vR(:,:,:,:,:,:,:,tR) = reshape(util,nk,na,na,na,ne,N,N);                          % last period indirect utility (vector!)
                                    end
                                % Initialize right-hand side of Bellman equation
                                    vmin=neg;
                                    l=0;
                                    
                                    % Loop over all k's in the capital grid to find the value,
                                    % which gives max of the right-hand side of Bellman equation
                                    
                                    while l < nk  	% assets tomorrow
                                        l=l+1;
                                        kap0 = kap(j); % current asset holdings
                                        kap1 = kap(l); % future asset holdings
                                        
                                        % Instantaneous utility
                                        cons1=(1+r0)*kap0-kap1;
                                        
                                        if cons1<=0
                                            util=neg;
                                        else
                                            util=(cons1^((1-sigma)*gamma))/(1-sigma);
                                        end
                                        
                                        % Right-hand side of Bellman equation
                                        v0 = util + beta*vR(l,i+1)+consump_grid_n(n_shock)+consump_grid_n(i_shock);
                                        
                                        % Store indirect utility and optimal saving 
                                        if v0>vmin
                                            vR(j,e,a_em,i)=v0;
                                            kapRopt(j,i)=l;
                                            vmin=v0;
                                        end        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Working households
    for i=tW:-1:1               % age
        for e=1:N               % productivity shock
            for j=1:nk          % assets today
                
                % Initialize right-hand side of Bellman equation
                vmin=neg;
                l=0;
                
                % Loop over all k's in the capital grid to find the value,
                % which gives max of the right-hand side of Bellman equation
                
                while l<nk  	% assets tomorrow
                    l = l+1;
                    
                    kap0 = kap(j); % current asset holdings
                    kap1 = kap(l); % future asset holdings
                    
                    % Optimal labor supply
                    lab = (gamma*(1-theta)*z(e)*eff(i)*w0-(1-gamma)*((1+r0)*kap0-kap1))/((1-theta)*w0*z(e)*eff(i));
                    
                    % Check feasibility of labor supply
                    if lab>1
                        lab=1;
                    elseif lab<0
                        lab=0;
                    end
                    
                    % Instantaneous utility
                    cons1=(1+r0)*kap0+(1-theta)*w0*z(e)*eff(i)*lab-kap1;
                    
                    if cons1<=0
                        util=neg;
                    else
                        util=(((cons1^gamma)*(1-lab)^(1-gamma))^(1-sigma))/(1-sigma);
                    end
                    
                    % Right-hand side of Bellman equation
                    
                    % Need to be careful in terms of matrices, when transiting from last period of
                    % work to retirement:
                    
                    if i==tW       % retired next period
                        v0=util + beta*vR(l,1);
                    else
                        v0=util + beta*(tran(e,1)*vW(1,l,i+1)+tran(e,2)*vW(2,l,i+1));
                    end
                    
                    % Store indirect utility, optimal saving and labor
                    if v0>vmin
                        vW(e,j,i)=v0;
                        kapWopt(e,j,i)=l;
                        labopt(e,j,i)=lab;
                        vmin=v0;
                    end
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aggregate capital stock and employment                                  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializations
kgen=zeros(J,1);        % Aggregate capital for each generation
gk=zeros(nk,J);         % Distribution of agents over capital for each cohort
gkR=zeros(nk,tR);       % Distribution of retirees over capital
   
% Distribution of newborns over capital and shocks (k=0).
gkW=zeros(N*nk,tW);     % Distribution of workers over capital and shocks for each working cohort
% We write it as a vector: [k1z1 k1z2 k2z1 k2z2 k3z1 k3z2 ...]
% This form will be very convenient when we compute the distribution below
gkW(1,1)=z_init(1)*mass(1);  
gkW(2,1)=z_init(2)*mass(1);

% Distribution of newborns over capital (k=0).
gk(1,1)=mass(1);    

% Aggregate labor supply by generation
labgen=zeros(tW,1);     % Aggregate labor for each generation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iterating over the distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workers
for d=1:tW             % iterations over cohort
    % Transition matrix
    % It is a square matrix. The ordering of rows and columns is the same:
    % [k1z1 k1z2 k2z1 k2z2 k3z1 k3z2 ...]
    clear trans
    trans=zeros(nk*N,nk*N);
    
   for ii = 1:nk; 
       for jj = 1:N; 
           trans((ii-1)*N+jj,(kapWopt(jj,ii,d)-1)*N+1:kapWopt(jj,ii,d)*N) = tran(jj,:)'; 
       end
   end
   
   % Aggregate labor suply by age
   gkW_resh=reshape(gkW(:,d),N,nk);
   
   labgen(d,1)=0;
   for ii = 1:nk;
       for jj = 1:N;
           labgen(d,1)=labgen(d,1)+z(jj)*labopt(jj,ii,d)*gkW_resh(jj,ii);
       end
   end
   
   % Aggregate capital by age
   if d<tW
       gkW(:,d+1)=trans'*gkW(:,d)/(1+n);
       gkW_resh=reshape(gkW(:,d+1),N,nk);
       gk(:,d+1)=sum(gkW_resh);
       kgen(d+1,1)=kap*gk(:,d+1);
       
   % Need to be careful when transiting from last period of work to retirement:
   elseif d==tW
       gkR_aux=trans'*gkW(:,d)/(1+n);
       gkR_aux_resh=reshape(gkR_aux,N,nk);
       gkR(:,1)=sum(gkR_aux_resh);
       
       gk(:,d+1)=gkR(:,1);
       kgen(d+1,1)=kap*gk(:,d+1);
   end
end
    
% Retirees
for d=1:tR-1   
    % Transition matrix
    % It is a square matrix. The ordering of rows and columns is the same:
    % [k1 k2 k3 ...]
    clear trans
    trans=zeros(nk,nk);
   
   for ii = 1:nk
       trans((ii-1)+1,(kapRopt(ii,d)-1)+1:kapRopt(ii,d)) = 1.0; 
   end

   gkR(:,d+1)=trans'*gkR(:,d)/(1+n);  
   kgen(d+1+tW,1)=kap*gkR(:,d+1); 
   gk(:,d+1+tW)=gkR(:,d+1);      
 
end  

k1=sum(kgen);
l1=sum(labgen);

%% Update the guess on capital and labor    
k0=0.95*k0+0.05*k1;
l0=0.9*l0+0.05*l1;

%% Display results
disp('  capital     labor   pension');
disp([k0, l0, ss0]);
disp('deviation-capital deviation-labor       ');
disp([abs(k1-k0),  abs(l1-l0)]);
end 

%% Display equilibrium results
disp('      k0         l0       w         r         ss0   ');
disp([k0, l0, w0, r0, ss0]);

% Prices
r0  = alpha*(k0^(alpha-1))*(l0^(1-alpha))-delta;
w0=(1-alpha)*(k0^(alpha))*(l0^(-alpha));

% Welfare
welfare=0;

for i=1:tW
    gkW_resh=reshape(gkW(:,i),N,nk);
    sum_aux=gkW_resh.*vW(:,:,i);
    welfare=welfare+sum(sum_aux(:));
end

for i=1:tR
    sum_aux=gkR.*vR;
    welfare=welfare+sum(sum_aux(:));
end

% Coefficient of variation    
sum_aux=(kap-k0).*(kap-k0);
sum_aux1=sum(gk');
variance=sum(sum_aux1.*sum_aux);
cv=sqrt(variance)/k0

%% Plots
% Value function for a retired agent
figure(1)
age=5;
plot1=plot(kap,vR(:,age));
set(plot1,'LineWidth', 1.5);
xlabel('capital','FontSize',14);
ylabel('value function','FontSize',14);
%title(['value function of a retired agent at age ', num2str(tW+age)])
shg
print('vf_retired','-dpng','-r300')

% Savings of a working agent
figure(2)
age=20;
plot1=plot(kap,kap(kapWopt(1,:,age)),'k-',kap,kap(kapWopt(2,:,age)),'r--',kap,kap,'-x');
set(plot1,'LineWidth', 1.0);
xlabel('private asset holdings','FontSize',14);
ylabel('saving','FontSize',14);
legend('high shock', 'low shock','45 degree')
%title(['saving of a working agent at age ', num2str(age)])
shg
print('saving_worker','-dpng','-r300')

% Savings profile
figure(3)
plot(1:J,kgen);
xlabel('age','FontSize',14);
ylabel('saving','FontSize',14);
title('savings profile of workers')
shg

% Labor supply of a working agent
figure(5)
age=20;
plot(kap,labopt(1,:,age),'k-',kap,labopt(2,:,age),'r--');
xlabel('capital','FontSize',14);
ylabel('saving','FontSize',14);
legend('high shock', 'low shock')
title(['labor supply of a working agent at age ', num2str(age)])
shg