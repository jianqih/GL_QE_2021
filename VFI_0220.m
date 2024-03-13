%% Guo_Leung_2021_QE 
clc
clear
close all

parameters;
% declare the ability grid

value_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); %# 10 dims
career_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); % 10 dims
capitalchoice_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK,SHOCK); %11 dims

value = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
career = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
capitalchoice = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK,SHOCK); % 11 dims
transfer = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
value_worker_noshock = ones(RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE,CAP); % 8 dims
value_non_incorp_noshock = ones(RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE,CAP); % 8 dims
value_incorp_noshock = ones(RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP); % 8 dims

schoolvalue = ones(ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolchoice = ones(ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolcapital = ones(ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolchoice2 = ones(ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolcapital2 = ones(ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
value_non_elite_noshock = ones(ABI,ABI,ABI,CAP,CAP); % 6 dims
value_elite_noshock = ones(ABI,ABI,ABI,CAP,CAP); % 6 dims

borrow_incorp = zeros(AGE,CAP);
borrow_non_incorp = zeros(AGE,CAP);
K_worker = zeros(AGE,CAP);
iter = 0;
err = 1;
tol_r = 0.000001;
maxrate = 0.20;
minrate = -0.20;

W = 1;

% Preparation
walue  = zeros(size(uM));
valueM = zeros(cS.nk*cS.ns,1,Ts.T-1);
Ind    = zeros(cS.nk*cS.ns,1,Ts.T-1);
vNext  = zeros(cS.ns, cS.nk, Ts.T-1);

while abs(err)>tol_r
    iter = iter + 1;
    R = (minrate+maxrate)*0.5; %
    fprintf("Iteration: %d\n",iter);
    for t = AGE:-1:RETIRE_AGE+1
        vNext(:,:,t)                = valueM_ts(:,:,t+1)';
        walue(:,:,t)                = uM(:,:,t) + cS.beta*kron(cS.P*vNext(:,:,t),ones(cS.nk,1));    
        [valueM(:,:,t), Ind(:,:,t)] = max(walue(:,:,t),[],2);
        valueM_ts(:,:,t)            = reshape(valueM(:,:,t),cS.nk,cS.ns);
        Ind_ts(:,:,t)               = reshape(Ind(:,:,t),cS.nk,cS.ns);  
    end
end
