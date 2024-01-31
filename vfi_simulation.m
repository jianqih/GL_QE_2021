%% Guo_Leung_2021_QE 
clc
clear
close all

parameters;
% declare the ability grid
abi_grid_em = randn(ABI,1);
abi_grid_ub = randn(ABI,1);
abi_grid_ib = randn(ABI,1);

for i=1:ABI
    abi_grid_em(i) = -RANGE_ABI * ABI_EM_STD + 2 * RANGE_ABI * ABI_EM_STD * i / (ABI - 1);
    abi_grid_ub(i) = -RANGE_ABI * ABI_UB_STD + 2 * RANGE_ABI * ABI_UB_STD * i / (ABI - 1);
    abi_grid_ib(i) = -RANGE_ABI * ABI_IB_STD + 2 * RANGE_ABI * ABI_IB_STD * i / (ABI - 1);
    fprintf("abi_grid_em[%d]=%f\n",i,abi_grid_em(i));
    fprintf("abi_grid_ub[%d]=%f\n",i,abi_grid_ub(i));
    fprintf("abi_grid_ib[%d]=%f\n",i,abi_grid_ib(i));
end

mid_abi_grid_em = ones([ABI-1,1]);
mid_abi_grid_ub = ones([ABI-1,1]);
mid_abi_grid_ib = ones([ABI-1,1]);
prob_abi_em = ones([ABI,1]);
prob_abi_ub = ones([ABI,1]);
prob_abi_ib = ones([ABI,1]);

for i=1:ABI-1
    mid_abi_grid_em(i) = (abi_grid_em(i) + abi_grid_em(i+1)) / 2;
    mid_abi_grid_ub(i) = (abi_grid_ub(i) + abi_grid_ub(i+1)) / 2;
    mid_abi_grid_ib(i) = (abi_grid_ib(i) + abi_grid_ib(i+1)) / 2;
    fprintf("mid_abi_grid_em[%d]=%f\n", i, mid_abi_grid_em(i));
    fprintf("mid_abi_grid_ub[%d]=%f\n", i, mid_abi_grid_ub(i));
    fprintf("mid_abi_grid_ib[%d]=%f\n", i, mid_abi_grid_ib(i));
end

prob_abi_em(1) = normcdf(mid_abi_grid_em(1), 0, ABI_EM_STD);
for i=2:ABI-1
    prob_abi_em(i) = normcdf(mid_abi_grid_em(i), 0, ABI_EM_STD) - normcdf(mid_abi_grid_em(i-1), 0, ABI_EM_STD);
end
prob_abi_em(ABI) = 1 - normcdf(mid_abi_grid_em(ABI-1), 0, ABI_EM_STD);

prob_abi_ub(1) = normcdf(mid_abi_grid_ub(1), 0, ABI_UB_STD);
for i=2:ABI-1
    prob_abi_ub(i) = normcdf(mid_abi_grid_ub(i), 0, ABI_UB_STD) - normcdf(mid_abi_grid_ub(i-1), 0, ABI_UB_STD);
end
prob_abi_ub(ABI) = 1 - normcdf(mid_abi_grid_ub(ABI-1), 0, ABI_UB_STD);

prob_abi_ib(1) = normcdf(mid_abi_grid_ib(1), 0, ABI_IB_STD);
for i=2:ABI-1
    prob_abi_ib(i) = normcdf(mid_abi_grid_ib(i), 0, ABI_IB_STD) - normcdf(mid_abi_grid_ib(i-1), 0, ABI_IB_STD);
end
prob_abi_ib(ABI) = 1 - normcdf(mid_abi_grid_ib(ABI-1), 0, ABI_IB_STD);

for i=1:ABI
    fprintf("prob_abi_em[%d]=%f\n", i, prob_abi_em(i));
    fprintf("prob_abi_ub[%d]=%f\n", i, prob_abi_ub(i));
    fprintf("prob_abi_ib[%d]=%f\n", i, prob_abi_ib(i));
end
% shock grid
% wage and investment
% consumption value of non-incorporated and incorporated
% consumption value of elite and non-elite college
shock_grid_w = randn(EDU,SHOCK);
shock_grid_n = randn([SHOCK,1]);
shock_grid_i = randn([SHOCK,1]);
consump_grid_n = randn([SHOCK,1]);
consump_grid_i = randn([SHOCK,1]);
consump_grid_c = randn([SHOCK,1]);
consump_grid_e = randn([SHOCK,1]);

for i=1:SHOCK
    shock_grid_w(1,i) = -RANGE_SHOCK * W_H_STD + 2 * RANGE_SHOCK * W_H_STD * i / (SHOCK - 1);
    shock_grid_w(2,i) = -RANGE_SHOCK * W_C_STD + 2 * RANGE_SHOCK * W_C_STD * i / (SHOCK - 1);
    shock_grid_w(3,i) = -RANGE_SHOCK * W_E_STD + 2 * RANGE_SHOCK * W_E_STD * i / (SHOCK - 1);
    shock_grid_n(i) = -RANGE_SHOCK * N_STD + 2 * RANGE_SHOCK * N_STD * i / (SHOCK - 1);
    shock_grid_i(i) = -RANGE_SHOCK * I_STD + 2 * RANGE_SHOCK * I_STD * i / (SHOCK - 1);
    consump_grid_n(i) = -RANGE_SHOCK * CONSUMP_N_STD + 2 * RANGE_SHOCK * CONSUMP_N_STD * i / (SHOCK - 1);
    consump_grid_i(i) = -RANGE_SHOCK * CONSUMP_I_STD + 2 * RANGE_SHOCK * CONSUMP_I_STD * i / (SHOCK - 1);
    consump_grid_c(i) = -RANGE_SHOCK * CONSUMP_C_STD + 2 * RANGE_SHOCK * CONSUMP_C_STD * i / (SHOCK - 1);
    consump_grid_e(i) = -RANGE_SHOCK * CONSUMP_E_STD + 2 * RANGE_SHOCK * CONSUMP_E_STD * i / (SHOCK - 1);
    fprintf("consump_grid_n[%d]=%f\n", i, consump_grid_n(i));
    fprintf("shock_grid_i[%d]=%f\n", i, shock_grid_i(i));
    
end
mid_shock_grid_w = ones(EDU,SHOCK - 1);
mid_shock_grid_n = ones([SHOCK - 1,1]);
mid_shock_grid_i = ones([SHOCK - 1,1]);
mid_consump_grid_n = ones([SHOCK - 1,1]);
mid_consump_grid_i = ones([SHOCK - 1,1]);
mid_consump_grid_c = ones([SHOCK - 1,1]);
mid_consump_grid_e = ones([SHOCK - 1,1]);
prob_shock_w = ones(EDU,SHOCK);
prob_shock_n = ones(SHOCK - 1,1);
prob_shock_i = ones(SHOCK - 1,1);
prob_consump_n = ones(SHOCK - 1,1);
prob_consump_i = ones(SHOCK - 1,1);
prob_consump_c = ones(SHOCK - 1,1);
prob_consump_e = ones(SHOCK - 1,1);

for i=1:SHOCK-1
    mid_shock_grid_w(1,i) = (shock_grid_w(1,i) + shock_grid_w(1,i+1)) / 2;
    mid_shock_grid_w(2,i) = (shock_grid_w(2,i) + shock_grid_w(2,i+1)) / 2;
    mid_shock_grid_w(3,i) = (shock_grid_w(3,i) + shock_grid_w(3,i)) / 2;
    mid_shock_grid_n(i) = (shock_grid_n(i) + shock_grid_n(i+1)) / 2;
    mid_shock_grid_i(i) = (shock_grid_i(i) + shock_grid_i(i+1)) / 2;
    mid_consump_grid_n(i) = (consump_grid_n(i) + consump_grid_n(i+1) ) / 2;
    mid_consump_grid_i(i) = (consump_grid_i(i) + consump_grid_i(i+1) ) / 2;
    mid_consump_grid_c(i) = (consump_grid_n(i) + consump_grid_n(i+1) ) / 2;
    mid_consump_grid_e(i) = (consump_grid_i(i) + consump_grid_i(i+1) ) / 2;
    fprintf("mid_shock_grid_w[%d]=%f\n", i, mid_shock_grid_w(i));
    fprintf("mid_shock_grid_i[%d]=%f\n", i, mid_shock_grid_i(i));
end

prob_shock_w(1,1) = normcdf(mid_shock_grid_w(1,1), 0, W_H_STD);
for i=2:SHOCK-1
    prob_shock_w(1,i) = normcdf(mid_shock_grid_w(1,i), 0, W_H_STD) - normcdf(mid_shock_grid_w(1,i-1), 0, W_H_STD);
end
prob_shock_w(1,SHOCK) = 1 - normcdf(mid_shock_grid_w(1,SHOCK-1), 0, W_H_STD);

prob_shock_w(2,1) = normcdf(mid_shock_grid_w(2,1), 0, W_C_STD);
for i=2:SHOCK-1
    prob_shock_w(2,i) = normcdf(mid_shock_grid_w(2,i), 0, W_C_STD) - normcdf(mid_shock_grid_w(2,i-1), 0, W_C_STD);
end
prob_shock_w(2,SHOCK) = 1 - normcdf(mid_shock_grid_w(1,SHOCK-1), 0, W_C_STD);

prob_shock_w(3,SHOCK) = normcdf(mid_shock_grid_w(3,1), 0, W_E_STD);
for i=2:SHOCK-1
    prob_shock_w(3,i) = normcdf(mid_shock_grid_w(3,i), 0, W_E_STD) - normcdf(mid_shock_grid_w(3,i-1), 0, W_E_STD);
end
prob_shock_w(3,SHOCK) = 1 - normcdf(mid_shock_grid_w(3,SHOCK-1), 0, W_E_STD);

prob_shock_n(1) = normcdf(mid_shock_grid_n(1), 0, N_STD);
for i=2:SHOCK-1
    prob_shock_n(i) = normcdf(mid_shock_grid_n(i), 0, N_STD) - normcdf(mid_shock_grid_n(i-1), 0, N_STD);
end
prob_shock_n(SHOCK) = 1 - normcdf(mid_shock_grid_n(SHOCK-1), 0, N_STD);

prob_shock_i(1) = normcdf(mid_shock_grid_i(1), 0, I_STD);
for i=2:SHOCK-1
    prob_shock_i(i) = normcdf(mid_shock_grid_i(i), 0, I_STD) - normcdf(mid_shock_grid_i(i-1), 0, I_STD);
end
prob_shock_i(SHOCK) = 1 - normcdf(mid_shock_grid_i(SHOCK-1), 0, I_STD);

prob_consump_n(1) = normcdf(mid_consump_grid_n(1), 0, CONSUMP_N_STD);
for i=2:SHOCK-1
    prob_consump_n(i) = normcdf(mid_consump_grid_n(i), 0, CONSUMP_N_STD) - normcdf(mid_consump_grid_n(i-1), 0, CONSUMP_N_STD);
end
prob_consump_n(SHOCK) = 1 - normcdf(mid_consump_grid_n(SHOCK-1), 0, CONSUMP_N_STD);

prob_consump_i(1) = normcdf(mid_consump_grid_i(1), 0, CONSUMP_I_STD);
for i=2:SHOCK-1
    prob_consump_i(i) = normcdf(mid_consump_grid_i(i), 0, CONSUMP_I_STD) - normcdf(mid_consump_grid_i(i-1), 0, CONSUMP_I_STD);
end
prob_consump_i(SHOCK) = 1 - normcdf(mid_consump_grid_i(SHOCK-1), 0, CONSUMP_I_STD);

prob_consump_c(1) = normcdf(mid_consump_grid_c(1), 0, CONSUMP_C_STD);
for i=2:SHOCK-1
    prob_consump_c(i) = normcdf(mid_consump_grid_c(i), 0, CONSUMP_C_STD) - normcdf(mid_consump_grid_c(i-1), 0, CONSUMP_C_STD);
end
prob_consump_c(SHOCK) = 1 - normcdf(mid_consump_grid_c(SHOCK-1), 0, CONSUMP_C_STD);

prob_consump_e(1) = normcdf(mid_consump_grid_e(1), 0, CONSUMP_E_STD);
for i=2:SHOCK-1
    prob_consump_e(i) = normcdf(mid_consump_grid_e(i), 0, CONSUMP_E_STD) - normcdf(mid_consump_grid_e(i-1), 0, CONSUMP_E_STD);
end
prob_consump_e(SHOCK) = 1 - normcdf(mid_consump_grid_e(SHOCK-1), 0, CONSUMP_E_STD);
for i=1:SHOCK
    fprintf("prob_shock_n[%d]=%f\n", i,prob_shock_n(i));
    fprintf("prob_shock_i[%d]=%f\n", i,prob_shock_i(i));
end   

% ability intergenerational transition probability
tranprob_abi_em = randn(ABI,ABI);
tranprob_abi_ub = randn(ABI,ABI);
tranprob_abi_ib = randn(ABI,ABI);
for i = 1:ABI-1
    % at the boundary:
    tranprob_abi_em(i,1) = normcdf(mid_abi_grid_em(1) - THETA_EM * mid_abi_grid_em(i), 0, ABI_EM_STD);
    tranprob_abi_em(i,ABI) = 1 - normcdf(mid_abi_grid_em(ABI-1) - THETA_EM * mid_abi_grid_em(i), 0, ABI_EM_STD);
    tranprob_abi_ub(i,1) = normcdf(mid_abi_grid_ub(1) - THETA_UB * mid_abi_grid_ub(i), 0, ABI_UB_STD);
    tranprob_abi_ub(i,ABI) = 1 - normcdf(mid_abi_grid_ub(ABI-1) - THETA_UB * mid_abi_grid_ub(i), 0, ABI_UB_STD);
    tranprob_abi_ib(i,1) = normcdf(mid_abi_grid_ib(1) - THETA_IB * mid_abi_grid_ib(i), 0, ABI_IB_STD);
    tranprob_abi_ib(i,ABI) = 1 - normcdf(mid_abi_grid_ib(ABI-1) - THETA_IB * mid_abi_grid_ib(i), 0, ABI_IB_STD);
    for j=2:ABI-1
        tranprob_abi_em(i,j) = normcdf(mid_abi_grid_em(j) - THETA_EM * mid_abi_grid_em(i), 0, ABI_EM_STD) - normcdf(mid_abi_grid_em(j-1) - THETA_EM * mid_abi_grid_em(i), 0, ABI_EM_STD);
        tranprob_abi_ub(i,j) = normcdf(mid_abi_grid_ub(j) - THETA_UB * mid_abi_grid_ub(i), 0, ABI_UB_STD) - normcdf(mid_abi_grid_ub(j-1) - THETA_UB * mid_abi_grid_ub(i), 0, ABI_UB_STD);
        tranprob_abi_ib(i,j) = normcdf(mid_abi_grid_ib(j) - THETA_IB * mid_abi_grid_ib(i), 0, ABI_IB_STD) - normcdf(mid_abi_grid_ib(j-1) - THETA_IB * mid_abi_grid_ib(i), 0, ABI_IB_STD);
    end
end

value_retire = ones(ITER,AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); %# 10 dims
career_retire = ones(ITER,AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); % 10 dims
capitalchoice_retire = ones(ITER,AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK,SHOCK); %11 dims

value = ones(ITER,RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
career = ones(ITER,RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
capitalchoice = ones(ITER,RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK,SHOCK); % 11 dims
transfer = ones(ITER,RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 10 dims
value_worker_noshock = ones(ITER,RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE,CAP); % 8 dims
value_non_incorp_noshock = ones(ITER,RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE,CAP); % 8 dims
value_incorp_noshock = ones(ITER,RETIRE_AGE,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP); % 8 dims

schoolvalue = ones(ITER,ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolchoice = ones(ITER,ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolcapital = ones(ITER,ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolchoice2 = ones(ITER,ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
schoolcapital2 = ones(ITER,ABI,ABI,ABI,CAP,CAP,SHOCK,SHOCK); % 8 dims
value_non_elite_noshock = ones(ITER,ABI,ABI,ABI,CAP,CAP); % 6 dims
value_elite_noshock = ones(ITER,ABI,ABI,ABI,CAP,CAP); % 6 dims
for iter=1:ITER
    fprintf("Iteration: %d\n",iter);
    for t = AGE:-1:RETIRE_AGE+1
        fprintf("Processing Retirement age\n");
        fprintf("Processing age = %d\n",age(t));
        for e=1:EDU
            for a_em = 1:ABI % ability of employee
                for a_ub = 1:ABI % abi of unincorporated 
                    for a_ib= 1:ABI
                        for exp_bus=1:EXP_BUS_TYPE_RETIRE
                            futurevalue_retire0 = zeros(CAP,1);
                            futurevalue_retire1 = zeros(CAP,1);
                            futurevalue_retire2 = zeros(CAP,1);
                            for k_next_idx = 1:CAP
                                if (t < AGE)
                                    for consump_n_next=1:SHOCK
                                        for consump_i_next=1:SHOCK
                                            futurevalue_retire0(k_next_idx) = futurevalue_retire0(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                            futurevalue_retire1(k_next_idx) = futurevalue_retire1(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,2,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                            futurevalue_retire2(k_next_idx) = futurevalue_retire2(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,3,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                        end
                                    end
                                end
                            end
                            for k = 1:CAP % index k 
                                % worker's decision
                                value_worker = -9999;
                                capitalchoice_worker = -1;
                                k_worker = capital(k) * (1 + R) + P * W * h_em(a_em, e) * exp((4 * ALPHA_1(e) + 22.7 * ALPHA_2(e))); % actural holding
                                if (t == AGE)
                                    value_worker = util(k_worker);
                                    capitalchoice_worker = 1; % polciy function: capital = 0 
                                else
                                    for k_next = 0:CAP_STEP:round(min(k_worker, capital(CAP))) % real capital k next period
                                        value_tmp = util(k_worker - k_next) + futurevalue_retire0((k_next - CAP_MIN)/CAP_STEP+1);
                                        % disp((k_next - CAP_MIN)/CAP_STEP+1);
                                        if (value_worker < value_tmp)
                                            value_worker = value_tmp;
                                            capitalchoice_worker = (k_next - CAP_MIN)/CAP_STEP+1;
                                        end
                                    end
                                end

                                % non-incorporated business owner decision
                                value_non_incorp = 0;
                                capitalchoice_non_incorp = ones(SHOCK,1);
                                % only able to run a business when there is no debt
                                if (capital(k) >= 0)
                                    for shock_n = 1:SHOCK
                                        value_non_incorp_tmp = -9999;
                                        % capitalchoice_non_incorp(shock_n) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                                        invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp .* NU_UB), 1/(NU_UB - 1));
                                        % liquidity constraint
                                        if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                                            invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB)- I_COST(2) + (capital(k) - invest_non_incorp) * (1 + R);
                                        if (t == AGE)
                                            value_non_incorp_tmp = util(k_non_incorp);
                                            capitalchoice_non_incorp(shock_n) = 1;
                                        else
                                            for k_next = 0:CAP_STEP:min(k_non_incorp, capital(CAP))
                                                value_tmp = util(k_non_incorp - k_next) + futurevalue_retire1((k_next - CAP_MIN)/CAP_STEP+1);
                                                if (value_non_incorp_tmp < value_tmp)
                                                    value_non_incorp_tmp = value_tmp;
                                                    capitalchoice_non_incorp(shock_n) = (k_next - CAP_MIN)/CAP_STEP+1;
                                                end
                                            end
                                        end
                                        value_non_incorp = value_non_incorp + value_non_incorp_tmp * prob_shock_n(shock_n);
                                    end
                                else
                                    value_non_incorp = -999;
                                end

                                % incorporated business owner decision
                                value_incorp = 0;
                                capitalchoice_incorp = ones(SHOCK,1);
                                % only able to run a business when there is no debt
                                if (capital(k) >= 0)
                                    for shock_i=1:SHOCK
                                        value_incorp_tmp = -9999;
                                        capitalchoice_incorp(shock_i) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                                        invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                                        % liquidity constraint
                                        if invest_incorp > (1 + LAMBDA_E) * capital(k)
                                            invest_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST(1) + (capital(k) - invest_incorp) * (1 + R);
                                        if (exp_bus == 2)
                                            k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) + (capital(k) - invest_incorp) * (1 + R);
                                        end
                                        if (t == AGE)
                                            value_incorp_tmp = util(k_incorp);
                                            capitalchoice_incorp(shock_i) = 1;
                                        else
                                            for k_next=0:CAP_STEP:min(k_incorp,capital(CAP))
                                                value_tmp = util(k_incorp - k_next) + futurevalue_retire2((k_next - CAP_MIN)/CAP_STEP+1);
                                                if (value_incorp_tmp < value_tmp)
                                                    value_incorp_tmp = value_tmp;
                                                    capitalchoice_incorp(shock_i) = (k_next - CAP_MIN)/CAP_STEP+1;
                                                end
                                            end
                                        end
                                        value_incorp = value_incorp + value_incorp_tmp * prob_shock_i(shock_i);
                                    end
                                else
                                    value_incorp = -999;
                                end

                                for consump_n = 1:SHOCK
                                    for consump_i = 1:SHOCK
                                        if (exp_bus > 1)
                                            if value_worker > value_non_incorp + consump_grid_n(consump_n) && value_worker > value_incorp + consump_grid_i(consump_i)
                                                value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                                                career_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                                                for shock_w = 1:SHOCK
                                                    capitalchoice_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
                                                end
                                            elseif value_non_incorp + consump_grid_n(consump_n) >= value_incorp + consump_grid_i(consump_i)
                                                value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_non_incorp + consump_grid_n(consump_n);
                                                career_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = NON_INCORP;
                                                for shock_n = 1:SHOCK
                                                    capitalchoice_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n) = capitalchoice_non_incorp(shock_n);
                                                end
                                            else
                                                value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_incorp + consump_grid_i(consump_i);
                                                career_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = INCORP;
                                                for shock_i = 1:SHOCK
                                                    capitalchoice_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i) = capitalchoice_incorp(shock_i);
                                                end
                                            end
                                        else
                                            value_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                                            career_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                                            for shock_w = 1:SHOCK
                                                capitalchoice_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
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
    end

    for t = RETIRE_AGE-1:-1:1
        fprintf("Processing Working age\n")
        fprintf("Processing age = %d\n",age(t))
    % pragma omp parallel for collapse(4)
        for e = 1:EDU
            for a_em= 1:ABI
                for a_ub = 1:ABI
                    for a_ib = 1:ABI
                        for exp_bus = 1:EXP_BUS_TYPE
                            futurevalue0 = zeros(CAP,1);
                            futurevalue1 = zeros(CAP,1);
                            futurevalue_retire0 = zeros(CAP,1);
                            futurevalue_retire1 = zeros(CAP,1);
                            futurevalue_retire2 = zeros(CAP,1);
                            for k_next_idx = 1:CAP
                                if (t == RETIRE_AGE)
                                    for consump_n_next = 1:SHOCK
                                        for consump_i_next= 1:SHOCK
                                            futurevalue_retire0(k_next_idx) = futurevalue_retire0(k_next_idx) + BETA * value_retire(iter,1,e,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                            futurevalue_retire1(k_next_idx) = futurevalue_retire0(k_next_idx) + BETA * value_retire(iter,1,e,a_em,a_ub,a_ib,2,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                            futurevalue_retire2(k_next_idx) = futurevalue_retire2(k_next_idx) + BETA * value_retire(iter,1,e,a_em,a_ub,a_ib,3,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                        end
                                    end
                                else
                                    for consump_n_next = 1:SHOCK
                                        for consump_i_next = 1:SHOCK
                                            futurevalue0(k_next_idx) = futurevalue0(k_next_idx) + BETA * value(iter,t+1,e,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                            futurevalue1(k_next_idx) = futurevalue0(k_next_idx) + BETA * value(iter,t+1,e,a_em,a_ub,a_ib,2,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                                        end
                                    end
                                end
                            end

                            for k=1:CAP
                                % worker's decision
                                value_worker = 0;
                                capitalchoice_worker = ones(SHOCK,1);
                                for shock_w = 1:SHOCK
                                    value_worker_tmp = -9999;
                                    capitalchoice_worker(shock_w) = -1;
                                    k_worker = capital(k) * (1 + R) + W * h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w)));
                                    for k_next = 0:CAP_STEP:round(min(k_worker,capital(CAP))) % capital cannot negative
                                        value_tmp = util(k_worker - k_next);
                                        if t == RETIRE_AGE
                                            value_tmp = value_tmp + futurevalue_retire0((k_next - CAP_MIN)/CAP_STEP+1);
                                        else 
                                            value_tmp = value_tmp + futurevalue0((k_next - CAP_MIN)/CAP_STEP+1);
                                        end
                                        if (value_worker_tmp < value_tmp)
                                            value_worker_tmp = value_tmp;
                                            capitalchoice_worker(shock_w) = (k_next - CAP_MIN)/CAP_STEP+1;
                                        end
                                    end
                                    value_worker = value_worker + value_worker_tmp * prob_shock_w(e,shock_w);
                                end
                                % non-incorporated business owner decision
                                value_non_incorp = 0;
                                capitalchoice_non_incorp = ones(SHOCK,1);
                                % only able to run a business when there is no debt
                                if (capital(k) >= 0)
                                    for shock_n = 1:SHOCK
                                        value_non_incorp_tmp = -9999;
                                        capitalchoice_non_incorp(shock_n) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                                        invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                                        % // liquidity constraint
                                        if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                                            invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) + (capital(k) - invest_non_incorp) * (1 + R);
                                        for k_next = 0:CAP_STEP:min(k_non_incorp,capital(CAP))
                                            value_tmp = util(k_non_incorp - k_next);
                                            if (t == RETIRE_AGE)
                                                value_tmp = value_tmp + futurevalue_retire1((k_next - CAP_MIN)/CAP_STEP+1);
                                            else
                                                value_tmp = value_tmp + futurevalue0((k_next - CAP_MIN)/CAP_STEP+1);
                                            end
                                            if (value_non_incorp_tmp < value_tmp)
                                                value_non_incorp_tmp = value_tmp;
                                                capitalchoice_non_incorp(shock_n) = (k_next - CAP_MIN)/CAP_STEP+1;
                                            end
                                        end
                                        value_non_incorp = value_non_incorp + value_non_incorp_tmp * prob_shock_n(shock_n);
                                        
                                    end
                                else
                                    value_non_incorp = -999;
                                end

                                % incorporated business owner decision
                                value_incorp = 0;
                                % double value_incorp_tmp;
                                capitalchoice_incorp = ones(SHOCK,1);
                                % only able to run a business when there is no debt
                                if (capital(k) >= 0)
                                    for shock_i =1:SHOCK
                                        value_incorp_tmp = -9999;
                                        capitalchoice_incorp(shock_i) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                                        invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                                        % liquidity constraint
                                        if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                                            invest_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST(1) + (capital(k) - invest_incorp) * (1 + R);
                                        if (exp_bus == 1)
                                            k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) + (capital(k) - invest_incorp) * (1 + R);
                                        end
                                        for k_next = 0:CAP_STEP:min(k_incorp,capital(CAP))
                                            value_tmp = util(k_incorp - k_next);
                                            if (t == RETIRE_AGE)
                                                value_tmp = value_tmp + futurevalue_retire2((k_next - CAP_MIN)/CAP_STEP+1);
                                            else
                                                value_tmp = value_tmp + futurevalue1((k_next - CAP_MIN)/CAP_STEP+1);
                                            end
                                            if (value_incorp_tmp < value_tmp)
                                                value_incorp_tmp = value_tmp;
                                                capitalchoice_incorp(shock_i) = (k_next - CAP_MIN)/CAP_STEP+1;
                                            end
                                        end
                                        value_incorp = value_incorp + value_incorp_tmp * prob_shock_i(shock_i);
                                    end
                                else
                                    value_incorp = -999;
                                end

                                
                                value_worker_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k) = value_worker;
                                value_non_incorp_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k) = value_non_incorp;
                                value_incorp_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k) = value_incorp;

                                for consump_n = 1:SHOCK
                                    for consump_i = 1:SHOCK
                                        if (value_worker > value_non_incorp + consump_grid_n(consump_n) && value_worker > value_incorp + consump_grid_i(consump_i))
                                            value(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                                            career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                                            for shock_w = 1:SHOCK
                                                capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k) = capitalchoice_worker(shock_w);
                                            end
                                        elseif (value_non_incorp + consump_grid_n(consump_n) >= value_incorp + consump_grid_i(consump_i))
                                            value(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_non_incorp + consump_grid_n(consump_n);
                                            career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = NON_INCORP;
                                            for shock_n = 1:SHOCK
                                                capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n) = capitalchoice_non_incorp(shock_n);
                                            end
                                        else
                                            value(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_incorp + consump_grid_i(consump_i);
                                            career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = INCORP;
                                            for shock_i = 1:SHOCK
                                                capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i) = capitalchoice_incorp(shock_i);
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
    end
    

    fprintf("start schooling decision\n");

    % schooling decision
    for a_em=1:ABI
        for a_ub=1:ABI
            for a_ib=1:ABI
                schoolfuture1 = zeros(CAP,1);
                schoolfuture2 = zeros(CAP,1);
                for k_next_idx=1:CAP
                    for consump_n_next=1:SHOCK
                        for consump_i_next=1:SHOCK
                            schoolfuture1(k_next_idx) = schoolfuture1(k_next_idx) + BETA * value(iter,2,2,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                            schoolfuture2(k_next_idx) = schoolfuture1(k_next_idx) + BETA * value(iter,2,3,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        end
                    end
                end

                for k = 1:CAP
                    for k_fam = 1:CAP
                        % high school edu = 1
                        value_highschool = 0;
                        for consump_n = 1:SHOCK
                            for consump_i = 1:SHOCK
                                value_highschool = value_highschool + value(iter,1,1,a_em,a_ub,a_ib,1,k,consump_n,consump_i) * prob_consump_n(consump_n) * prob_consump_i(consump_i);
                            end
                        end

                        % non-elite college edu = 2
                        value_college = -9999;
                        capitalchoice_college = -1;
                        k_college = (1 + R) * (capital(k) - COLLEGE_TUITION + financial_aid_nonelite(k_fam, a_em));
                        for k_next = CAP_MIN:CAP_STEP:min(k_college,capital(CAP)) % allow debt at college
                            value_tmp = util(k_college - k_next) + CONSUMPTION_COLLEGE + schoolfuture1((k_next - CAP_MIN)/CAP_STEP+1);
                            if (value_college < value_tmp)
                                value_college = value_tmp;
                                capitalchoice_college = round((k_next - CAP_MIN)/CAP_STEP+1);
                            end
                        end

                        % elite college = 3
                        value_elitecollege = -9999;
                        capitalchoice_elitecollege = -1;
                        k_elitecollege = (1 + R) * (capital(k) - ELITE_TUITION + financial_aid_elite(k_fam, a_em));
                        for k_next = CAP_MIN:CAP_STEP:min(k_elitecollege,capital(CAP))
                            value_tmp = util(k_elitecollege - k_next) + CONSUMPTION_ELITE + schoolfuture2((k_next - CAP_MIN)/CAP_STEP+1);
                            if (value_elitecollege < value_tmp)
                                value_elitecollege = value_tmp;
                                capitalchoice_elitecollege = (k_next - CAP_MIN)/CAP_STEP+1;
                            end
                        end

                        value_non_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam) = value_college;
                        value_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam) = value_elitecollege;

                        for consump_c = 1:SHOCK
                            for consump_e = 1:SHOCK
                                if (value_highschool > value_college + consump_grid_c(consump_c) && value_highschool > value_elitecollege + consump_grid_e(consump_e))
                                    schoolchoice(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = NO_COLLEGE;
                                    schoolvalue(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = value_highschool;
                                    schoolcapital(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = -1;
                                elseif (value_college + consump_grid_c(consump_c) > value_elitecollege + consump_grid_e(consump_e))
                                    schoolchoice(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = COLLEGE;
                                    schoolvalue(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = value_college + consump_grid_c(consump_c);
                                    schoolcapital(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = capitalchoice_college;
                                else
                                    schoolchoice(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = ELITE;
						            % check if the student is admitted by elite college
                                    if a_em >= 5
							            sat_g = 3;
                                    elseif a_em >= 2
						 	            sat_g = 2;
                                    else 
                                        sat_g = 1;
                                    end
                                    schoolvalue(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = ELITE_ADMIT(sat_g) * value_elitecollege + (1-ELITE_ADMIT(sat_g)) * value_college + consump_grid_e(consump_e);
                                    schoolcapital(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = capitalchoice_elitecollege;
                                    if (value_highschool > value_college + consump_grid_c(consump_c))
                                        schoolchoice2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = NO_COLLEGE;
                                        schoolcapital2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = -1;
                                    else
                                        schoolchoice2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = COLLEGE;
                                        schoolcapital2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = capitalchoice_college;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    % return;
    %% Simulation %%
    fprintf("start simulation\n");
    dim = 3;
    rng(1);
    mean = [0;0;0];
    covar = diag([ABI_EM_STD^2, ABI_UB_STD^2, ABI_IB_STD^2]);
    % Cholesky decomposition
    try
        cholSolver = chol(covar, 'lower');
        normTransform = cholSolver;
    catch
        % Eigen solver
        [eigenVectors, eigenValues] = eig(covar);
        sqrtEigenValues = sqrt(diag(eigenValues));
        normTransform = eigenVectors * diag(sqrtEigenValues);
    end
    
    % Generate random samples
    randN = randn(dim, num_id);
    samples = bsxfun(@plus, normTransform * randN, mean);
    
    % simulation array
    abi_em_sim = ones([num_id,1]);
    abi_ub_sim = ones([num_id,1]);
    abi_ib_sim = ones([num_id,1]);
    abi_em_org_sim = ones([num_id,1]);
    abi_ub_org_sim = ones([num_id,1]);
    abi_ib_org_sim = ones([num_id,1]);
    edu_sim = ones([num_id,1]);
    k_idx_sim = ones([num_id,AGE+1]); % individual capital simulation
    k_fam_idx_sim = ones([num_id,1]);
    career_sim = ones([num_id,AGE]);
    income_sim = ones([num_id,AGE]);
    exp_bus_sim = ones([num_id,AGE]);
    shock_w_sim = ones([num_id,AGE]);
    shock_n_sim = ones([num_id,AGE]);
    shock_i_sim = ones([num_id,AGE]);
    shock_s_sim = ones([num_id,AGE]);
    consump_n_sim = ones([num_id,AGE]);
    consump_i_sim = ones([num_id,AGE]);
    consump_c_sim = ones([num_id,AGE]);
    consump_e_sim = ones([num_id,AGE]);
    
    
    for id = 1:num_id
        % initial ability
        abi_em_rand = samples(1, id);
        abi_em_org_sim(id) = abi_em_rand;
        abi_em_min = abi_grid_em(1);
        abi_em_step = 2 * RANGE_ABI * ABI_EM_STD / (ABI - 1);
        abi_em_sim(id) = round((abi_em_rand - abi_em_min)/abi_em_step+1);
        if (abi_em_sim(id) > ABI)
            abi_em_sim(id) = ABI;
        elseif (abi_em_sim(id) < 1)
            abi_em_sim(id) = 1;
        end
        abi_ub_rand = samples(1, id);
        abi_ub_org_sim(id) = abi_ub_rand;
        abi_ub_min = abi_grid_ub(1);
        abi_ub_step = 2 * RANGE_ABI * ABI_UB_STD / (ABI - 1);
        abi_ub_sim(id) = round((abi_ub_rand - abi_ub_min)/abi_ub_step+1);
        if (abi_ub_sim(id) > ABI)
            abi_ub_sim(id) = ABI;
        elseif (abi_ub_sim(id) < 1)
            abi_ub_sim(id) = 1;
        end
        abi_ib_rand = samples(2, id);
        abi_ib_org_sim(id) = abi_ib_rand;
        abi_ib_min = abi_grid_ib(1);
        abi_ib_step = 2 * RANGE_ABI * ABI_IB_STD / (ABI - 1);
        abi_ib_sim(id) = round((abi_ib_rand - abi_ib_min)/abi_ib_step+1);
        if (abi_ib_sim(id) > ABI) % censored
            abi_ib_sim(id) = ABI;
        elseif (abi_ib_sim(id) < 1)
            abi_ib_sim(id) = 1;
        end
        % initial wealth
        k_random_number = (rand() / (RAND_MAX)) * 100; %//draw a random number between 0 and 100
        cdf = K_CON_DIST(abi_em_sim(id),1);
        group_id = 1;
        while (cdf < k_random_number)
            group_id = group_id + 1;
            cdf = cdf + K_CON_DIST(abi_em_sim(id), group_id);
        end
        k_idx_sim(id,1) = group_mapping(group_id,num_id,CAP); % initial wealth 
        % initial family wealth
        k_fam_rand = r8_normal(0, K_FAM_STD);
        if (capital(k_idx_sim(id,1)) > 0)
            k_fam = exp(log(capital(k_idx_sim(id,1))) * COEF_K + COEF_ABI(abi_em_sim(id)) + k_fam_rand);
        else
            k_fam = exp(COEF_ABI(abi_em_sim(id)) + k_fam_rand);
        end
        k_fam_idx_sim(id) = round((k_fam - CAP_MIN)/CAP_STEP+1);
        % if (id < 500)
        %    std::cout << k_idx_sim[id][0] << ", " << k_fam_rand << ", " << k_fam << ", " << k_fam_idx_sim[id] << std::endl;
        if (k_fam_idx_sim(id) < 1) % censored
            k_fam_idx_sim(id) = 1;
        elseif (k_fam_idx_sim(id) > CAP)
            k_fam_idx_sim(id) = CAP;
        end
        % if (id < 10)
        %    fout << abi_em_rand << ", " << abi_ub_rand << ", " << k_random_number << ", " << k_fam_rand << std::endl;
    end
    
    % Initial ability distribution
    sum_abi_em = zeros(ABI,1);
    sum_abi_ub = zeros(ABI,1);
    sum_abi_ib = zeros(ABI,1);
    
    for id=1:num_id
        a_em = abi_em_sim(id);
        a_ub = abi_ub_sim(id);
        a_ib = abi_ib_sim(id);
        sum_abi_em(a_em) = sum_abi_em(a_em) + 1;
        sum_abi_ub(a_ub) = sum_abi_ub(a_ub) + 1;
        sum_abi_ib(a_ib) = sum_abi_ib(a_ib) + 1;
    end
    
    frac_abi_em = sum_abi_em / num_id;
    frac_abi_ub = sum_abi_ub / num_id;
    frac_abi_ib = sum_abi_ib / num_id;
    
    % Initial wealth distribution
    sum_wealth = zeros(CAP,1);
    for id=1:num_id
        k = k_idx_sim(id,1);
        sum_wealth(k) = sum_wealth(k) + 1;
    end
    frac_wealth = sum_wealth / num_id;
    % Initial family wealth distribution
    sum_wealth_fam = zeros(CAP);
    for id = 1:num_id
        k = k_fam_idx_sim(id);
        sum_wealth_fam(k) = sum_wealth_fam(k) + 1;
    end
    frac_wealth_fam = sum_wealth_fam / num_id;
    % simulate shocks
    for id=1:num_id
        for t =1:RETIRE_AGE-1
            shock_w_rand = r8_normal(0, W_H_STD);
            shock_w_min = shock_grid_w(1,1);
            shock_w_step = 2 * RANGE_SHOCK * W_H_STD / (SHOCK - 1);
            shock_w_sim(id,t) = round((shock_w_rand - shock_w_min)/shock_w_step+1);
            if (shock_w_sim(id,t) > SHOCK)
                shock_w_sim(id,t) = SHOCK;
            elseif (shock_w_sim(id,t) < 1)
                shock_w_sim(id,t) = 1;
            end
            shock_n_rand = r8_normal(0, N_STD);
            shock_n_min = shock_grid_n(1);
            shock_n_step = 2 * RANGE_SHOCK * N_STD / (SHOCK - 1);
            shock_n_sim(id,t) = round((shock_n_rand - shock_n_min)/shock_n_step+1);
            if (shock_n_sim(id,t) > SHOCK)
                shock_n_sim(id,t) = SHOCK;
            elseif (shock_n_sim(id,t) < 1)
                shock_n_sim(id,t) = 1;
            end
            shock_i_rand = r8_normal(0, I_STD);
            shock_i_min = shock_grid_i(1);
            shock_i_step = 2 * RANGE_SHOCK * I_STD / (SHOCK - 1);
            shock_i_sim(id,t) = round((shock_i_rand - shock_i_min)/shock_i_step+1);
            if (shock_i_sim(id,t) > SHOCK)
                shock_i_sim(id,t) = SHOCK;
            elseif (shock_i_sim(id,t) < 1)
                shock_i_sim(id,t) = 1;
            end
            consump_n_rand = r8_normal(0, CONSUMP_N_STD);
            consump_n_min = consump_grid_n(1);
            consump_n_step = 2 * RANGE_SHOCK * CONSUMP_N_STD / (SHOCK-1);
            consump_n_sim(id,t) = round((consump_n_rand - consump_n_min)/consump_n_step+1);
            if (consump_n_sim(id,t) > SHOCK)
                consump_n_sim(id,t) = SHOCK;
            elseif (consump_n_sim(id,t) < 1)
                consump_n_sim(id,t) = 1;
            end
            consump_i_rand = r8_normal(0, CONSUMP_I_STD);
            consump_i_min = consump_grid_i(1);
            consump_i_step = 2 * RANGE_SHOCK * CONSUMP_I_STD / (SHOCK-1);
            consump_i_sim(id,t) = round((consump_i_rand - consump_i_min)/consump_i_step+1);
            if (consump_i_sim(id,t) > SHOCK)
                consump_i_sim(id,t) = SHOCK;
            elseif (consump_i_sim(id,t) < 1)
                consump_i_sim(id,t) = 1;
            end
        end
        consump_c_rand = r8_normal(0, CONSUMP_C_STD);
        consump_c_min = consump_grid_c(1);
        consump_c_step = 2 * RANGE_SHOCK * CONSUMP_C_STD / (SHOCK-1);
        consump_c_sim(id) = round((consump_c_rand - consump_c_min)/consump_c_step+1);
        if (consump_c_sim(id) > SHOCK)
            consump_c_sim(id) = SHOCK;
        elseif (consump_c_sim(id) < 1)
            consump_c_sim(id) = 1;
        end
        consump_e_rand = r8_normal(0, CONSUMP_E_STD);
        consump_e_min = consump_grid_e(1);
        consump_e_step = 2 * RANGE_SHOCK * CONSUMP_E_STD / (SHOCK-1);
        consump_e_sim(id) = round((consump_e_rand - consump_e_min)/consump_e_step+1);
        if (consump_e_sim(id) > SHOCK)
            consump_e_sim(id) = SHOCK;
        elseif (consump_e_sim(id) < 1)
            consump_e_sim(id) = 1;
        end
    
    end
    
    for id = 1:num_id
        for t = RETIRE_AGE:AGE-1
            shock_n_rand = r8_normal(0, N_STD);
            shock_n_min = shock_grid_n(1);
            shock_n_step = 2 * RANGE_SHOCK * N_STD / (SHOCK-1);
            shock_n_sim(id,t) = round((shock_n_rand - shock_n_min)/shock_n_step+1);
            if (shock_n_sim(id,t) > SHOCK)
                shock_n_sim(id,t) = SHOCK;
            elseif (shock_n_sim(id,t) < 1)
                shock_n_sim(id,t) = 1;
            end
            shock_i_rand = r8_normal(0, I_STD);
            shock_i_min = shock_grid_i(1);
            shock_i_step = 2 * RANGE_SHOCK * I_STD / (SHOCK - 1);
            shock_i_sim(id,t) = round((shock_i_rand - shock_i_min)/shock_i_step+1);
            if (shock_i_sim(id,t) > SHOCK)
                shock_i_sim(id,t) = SHOCK;
            elseif (shock_i_sim(id,t) < 1)
                shock_i_sim(id,t) = 1;
            end
            consump_n_rand = r8_normal(0, CONSUMP_N_STD);
            consump_n_min = consump_grid_n(1);
            consump_n_step = 2 * RANGE_SHOCK * CONSUMP_N_STD / (SHOCK - 1);
            consump_n_sim(id,t) = round((consump_n_rand - consump_n_min)/consump_n_step+1);
            if (consump_n_sim(id,t) > SHOCK)
                consump_n_sim(id,t) = SHOCK;
            elseif (consump_n_sim(id,t) < 1)
                consump_n_sim(id,t) = 1;
            end
            consump_i_rand = r8_normal(0, CONSUMP_I_STD);
            consump_i_min = consump_grid_i(1);
            consump_i_step = 2 * RANGE_SHOCK * CONSUMP_I_STD / (SHOCK - 1);
            consump_i_sim(id,t) = round((consump_i_rand - consump_i_min)/consump_i_step+1);
            if (consump_i_sim(id,t) > SHOCK)
                consump_i_sim(id,t) = SHOCK;
            elseif (consump_i_sim(id,t) < 1)
                consump_i_sim(id,t) = 1;
            end
            shock_s_sim(id,t) = (randn()/ RAND_MAX);
        end
    end
    
    laborw_sum = 0;
    capital_sum = 0;
    labor_sum = 0;
    capitalr_sum = 0;
    capitalw_sum = 0;
    % simulate first period: education decision
    fprintf("simulation first period\n");
    for id = 1:num_id
        iter = 1;
        t = 1;
        a_em = abi_em_sim(id);
        a_ub = abi_ub_sim(id);
        a_ib = abi_ib_sim(id);
        k = k_idx_sim(id,1);
        k_fam = k_fam_idx_sim(id);
        consump_c = consump_c_sim(id);
        consump_e = consump_e_sim(id);
        edu_random_number = (randn() / (RAND_MAX)); 
    
        edu_sim(id) = schoolchoice(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e);
        % check if the student is admitted by elite college
        flag = 0;
        if (edu_sim(id) == 3)
            sat_group=1;
	        if a_em >= 5
	            sat_group = 3;
            elseif a_em >= 2
	            sat_group = 2;
            end
	        if edu_random_number > ELITE_ADMIT(sat_group)
		        edu_sim(id) = schoolchoice2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e);
		        flag = 1; % admitted successfully.
            end  
        end
        
        % high school - work
        if (edu_sim(id) == 1)
            e = edu_sim(id); % e = 1
            exp_bus = 1; % without any business exp.
            shock_w = shock_w_sim(id,t);
            shock_n = shock_n_sim(id,t);
            shock_i = shock_i_sim(id,t);
            consump_n = consump_n_sim(id,t);
            consump_i = consump_i_sim(id,t);
    
            career_sim(id,t) = career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i);
            if (career_sim(id,t) == 1)
                % worker
                income_sim(id,t) = W * h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
                exp_bus_sim(id,t) = 1;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                labor_sum = labor_sum + h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
                laborw_sum = laborw_sum + income_sim(id,t);
                capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                capital_sum = capital_sum + capital(k_idx_sim(id,t));
                capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
            elseif (career_sim(id,1) == 2)
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                % liquidity constraint
                if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim(id,1) = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
                exp_bus_sim(id,1) = 2;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
                capital_sum = capital_sum + (invest_non_incorp - capital(k_idx_sim(id,t)));
                capitalr_sum = capitalr_sum + R * (invest_non_incorp - capital(k_idx_sim(id,t)));
            else
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                % liquidity constraint
                if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim(id,1) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST(1) - invest_incorp * (1 + R);
                exp_bus_sim(id,2) = 3;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
                capital_sum = capital_sum + (invest_incorp - capital(k_idx_sim(id,t)));
                capitalr_sum = capitalr_sum + R * (invest_incorp - capital(k_idx_sim(id,t)));
            end
        else
            % college
            career_sim(id,1) = -1;
            income_sim(id,1) = 1;
            exp_bus_sim(id,1) = 1;
            k_idx_sim(id,t+1) = schoolcapital(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e);
            if (flag == 1)
	            k_idx_sim(id,t+1) = schoolcapital2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e);
            end
        capital_sum = capital_sum + capital(k_idx_sim(id,t));
        capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
        end
    end
    
    % simulate period 2 - RETIRE_AGE
    fprintf("simulation second period to retire age\n");
    
    for id = 1:num_id
        for t = 2:RETIRE_AGE
            iter = 1;
            a_em = abi_em_sim(id);
            a_ub = abi_ub_sim(id);
            a_ib = abi_ib_sim(id);
            e = edu_sim(id);
            exp_bus = exp_bus_sim(id,t-1);
            k = k_idx_sim(id,t);
            consump_n = consump_n_sim(id,t);
            consump_i = consump_i_sim(id,t);
            shock_w = shock_w_sim(id,t);
            shock_n = shock_n_sim(id,t);
            shock_i = shock_i_sim(id,t);
    
            career_sim(id,t) = career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i);
            if (career_sim(id,t) == 1) % worker
                income_sim(id,t) = W * h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w)));
                exp_bus_sim(id,t) = 1;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                labor_sum = labor_sum + h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
                laborw_sum = laborw_sum + income_sim(id,t);
                capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                capital_sum = capital_sum + capital(k_idx_sim(id,t));
                capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
            elseif (career_sim(id,t) == 2)
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                % liquidity constraint
                if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
                exp_bus_sim(id,t) = 1;
                if (t == RETIRE_AGE)
                    exp_bus_sim(id,t) = 2;
                end
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
                capital_sum = capital_sum + (invest_non_incorp - capital(k_idx_sim(id,t)));
                capitalr_sum = capitalr_sum + R * (invest_non_incorp - capital(k_idx_sim(id,t)));
            else
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                % liquidity constraint
                if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST(1) - invest_incorp * (1 + R);
                if (exp_bus == 1)
                    income_sim(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
                end
                exp_bus_sim(id,t) = 1;
                if (t == RETIRE_AGE)
                    exp_bus_sim(id,t) = 2;
                end
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
                capital_sum = capital_sum + (invest_incorp - capital(k_idx_sim(id,t)));
                capitalr_sum = capitalr_sum + R * (invest_incorp - capital(k_idx_sim(id,t)));
            end
            
        end
    end
    
    % simulate RETIRE_AGE to Max_age
    fprintf("simulation retire age to max age\n");
    for id = 1:num_id
        for t = RETIRE_AGE+1:AGE
            iter = 1;
            a_em = abi_em_sim(id);
            a_ub = abi_ub_sim(id);
            a_ib = abi_ib_sim(id);
            e = edu_sim(id);
            exp_bus = exp_bus_sim(id,t-1);
            k = k_idx_sim(id,t);
            consump_n = consump_n_sim(id,t);
            consump_i = consump_i_sim(id,t);
            shock_w = shock_w_sim(id,t);
            shock_n = shock_n_sim(id,t);
            shock_i = shock_i_sim(id,t);
    
            if (shock_s_sim(id,t) < ZETA(t-RETIRE_AGE+1)) % survive or not
                career_sim(id,t) = career_retire(iter,t-RETIRE_AGE+1,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i);
                if (career_sim(id,t) == 1)
                    income_sim(id,t) = P * W * h_em(a_em, e) * exp(4 * ALPHA_1(e) + 22.7 * ALPHA_2(e));
                    exp_bus_sim(id,t) = 1;
                    k_idx_sim(id,t+1) = capitalchoice_retire(iter,t-RETIRE_AGE+1,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                    capital_sum = capital_sum + capital(k_idx_sim(id,t));
                    capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                    capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
                elseif (career_sim(id,t) == 2)
                    h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                    i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                    invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                    % liquidity constraint
                    if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                        invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    income_sim(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
                    exp_bus_sim(id,t) = 2;
                    k_idx_sim(id,t+1) = capitalchoice_retire(iter,t-RETIRE_AGE+1,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
                    capital_sum = capital_sum + (invest_non_incorp - capital(k_idx_sim(id,t)));
                    capitalr_sum = capitalr_sum + R * (invest_non_incorp - capital(k_idx_sim(id,t)));
                else
                    h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                    i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                    invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                    % liquidity constraint
                    if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                        invest_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    income_sim(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST(1) - invest_incorp * (1 + R);
                    if (exp_bus == 3)
                        income_sim(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
                    end
                    exp_bus_sim(id,t) = 3;
                    k_idx_sim(id,t+1) = capitalchoice_retire(iter,t-RETIRE_AGE+1,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
                    capital_sum = capital_sum + (invest_incorp - capital(k_idx_sim(id,t)));
                    capitalr_sum =  capitalr_sum + R * (invest_incorp - capital(k_idx_sim(id,t)));
                 end
            else
                career_sim(id,t) = -1;
            end
        end
    end
    % **********************************Calculate production function******************************************************************** 
    fprintf("production function\n");
    prod_alpha = capitalr_sum/(capitalr_sum + laborw_sum);
    prod_a = (R / prod_alpha) * power((labor_sum/capital_sum), prod_alpha - 1);
    fprintf("prod_alpha = %f\n prod_a = %f\n capitalr_sum = %f\n capital_sum = %f\n laborw_sum = %f\n labor_sum = %f\n", prod_alpha,prod_a,capitalr_sum, capital_sum,laborw_sum,labor_sum);
    W = PROD_EM.*alpha.*(labor_sum./(capitalw_sum-capital_sum)).^(1-alpha);
    R = PROD_EM.*(1-alpha).*(labor_sum./(capitalw_sum-capital_sum)).^(-alpha);
    fprintf("wage:%f\n interest rate:%f\n",W,R)
end
return;



% ************************************** Moment Calculation **************************************//
fprintf("finish simulation, start moment calculation\n");
% Proportion by edu group
fprintf("education choice\n");
sum_edu = zeros(EDU,1);
for id = 1:num_id
    e = edu_sim(id);
    sum_edu(e) = sum_edu(e) + 1;
end
frac_edu = sum_edu / num_id;
sum_abi_wealth = zeros(ABI,6);
sum_edu_abi_wealth = zeros(ABI,6,EDU);

for id = 1:num_id
    k_group = capital(k_idx_sim(id,1))/20000;
    if k_group < 1
        k_group = 1;
    elseif (k_group > 5)
        k_group = 5;
    end
    e = edu_sim(id);
    sum_abi_wealth(abi_em_sim(id),k_group) = sum_abi_wealth(abi_em_sim(id),k_group) + 1;
    sum_edu_abi_wealth(abi_em_sim(id),k_group,e) = sum_edu_abi_wealth(abi_em_sim(id),k_group,e) + 1;
end

% Ability by edu group
fprintf("ability by education group\n");
sum_abi_edu_em = zeros(EDU,1);
frac_abi_edu_em = zeros(EDU,1);
sum_abi_edu_ub = zeros(EDU,1);
frac_abi_edu_ub = zeros(EDU,1);
sum_abi_edu_ib = zeros(EDU,1);
frac_abi_edu_ib = zeros(EDU,1);
for id = 1:num_id
    e = edu_sim(id);
    sum_abi_edu_em(e) = sum_abi_edu_em(e) + abi_grid_em(abi_em_sim(id));
    sum_abi_edu_ub(e) = sum_abi_edu_ub(e) + abi_grid_ub(abi_ub_sim(id));
    sum_abi_edu_ib(e) = sum_abi_edu_ib(e) + abi_grid_ib(abi_ib_sim(id));
end
for e = 1:EDU
    frac_abi_edu_em(e) = sum_abi_edu_em(e) ./ sum_edu(e);
    frac_abi_edu_ub(e) = sum_abi_edu_ub(e) ./ sum_edu(e);
    frac_abi_edu_ib(e) = sum_abi_edu_ib(e) ./ sum_edu(e);
end

% Ability by edu group and career type
fprintf("ability by education group and career type\n");
sum_abi_educar_em = zeros(EDU,CAREER_CHOICES);
frac_abi_educar_em = zeros(EDU,CAREER_CHOICES);
sum_abi_educar_ub = zeros(EDU,CAREER_CHOICES);
frac_abi_educar_ub = zeros(EDU,CAREER_CHOICES);
sum_abi_educar_ib = zeros(EDU,CAREER_CHOICES);
frac_abi_educar_ib = zeros(EDU,CAREER_CHOICES);
count_ec = zeros(EDU,CAREER_CHOICES);

for id = 1:num_id
    for t = 2:8
        e = edu_sim(id);
        car = career_sim(id,t);
        sum_abi_educar_em(e,car) = sum_abi_educar_em(e,car) + abi_grid_em(abi_em_sim(id));
        sum_abi_educar_ub(e,car) = sum_abi_educar_ub(e,car) + abi_grid_ub(abi_ub_sim(id));
        sum_abi_educar_ib(e,car) = sum_abi_educar_ib(e,car) + abi_grid_ib(abi_ib_sim(id));
        count_ec(e,car) = count_ec(e,car) + 1;
    end
end
for e = 1:EDU
    for car = 1:CAREER_CHOICES
        frac_abi_educar_em(e,car) = sum_abi_educar_em(e,car)  / count_ec(e,car);
        frac_abi_educar_ub(e,car) = sum_abi_educar_ub(e,car)  / count_ec(e,car);
        frac_abi_educar_ib(e,car) = sum_abi_educar_ib(e,car)  / count_ec(e,car);
    end
end

% Proportion by career type
disp("career choice");
sum_career = zeros(CAREER_CHOICES,1);
count_career = 0;

for id = 1:num_id
    for t=1:RETIRE_AGE - 1
        car = career_sim(id,t);
        if (car >=1)
            sum_career(car) = sum_career(car) + 1;
            count_career = count_career + 1;
        end
    end
end
frac_career = sum_career./count_career;
% Career choice by age and edu group
sum_career_ageedu = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);
frac_career_ageedu = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);

for id = 1:num_id
    for t=1:RETIRE_AGE
        e = edu_sim(id);
        car = career_sim(id,t);
        if (car >=1)
            sum_career_ageedu(e,t,car) = sum_career_ageedu(e,t,car) + 1;
        end
    end
end

for e = 1:EDU
    for t = 2:RETIRE_AGE
        for car = 1:CAREER_CHOICES
            frac_career_ageedu(e,t,car) = sum_career_ageedu(e,t,car)/sum_edu(e);
        end
    end
end
% Career choice by edu group
fprintf("career choice by edu\n");

sum_career_edu = zeros(EDU,CAREER_CHOICES);
for id = 1:num_id
    for t = 2:RETIRE_AGE - 1
        e = edu_sim(id);
        car = career_sim(id,t);
        if (car >=1)
            sum_career_edu(e,car) = sum_career_edu(e,car) + 1;
        end
    end
end

frac_career_edu = sum_career_edu./sum_edu./(RETIRE_AGE - 2);

% Wage/Income by edu
fprintf("mean and sd of income by edu\n");
sum_income_edu = zeros(EDU,1);
sum_incomesq_edu = zeros(EDU,1);
for id = 1:num_id
    for t = 2:8
        e = edu_sim(id);
        sum_income_edu(e) = sum_income_edu(e) + income_sim(id,t);
    end
end
% avg_income_edu = sum_income_edu/sum_
avg_income_edu = sum_income_edu./sum_edu./(8-1);
for id = 1:num_id
    for t = 2:8
        e = edu_sim(id);
        sum_incomesq_edu(e) = sum_incomesq_edu(e) + power(income_sim(id,t) - avg_income_edu(e), 2);
    end
end
sd_income_edu = sqrt(sum_incomesq_edu./(sum_edu*(8-1)-1));

% Wage/Income by age
fprintf("average income by age\n");
sum_income_age = zeros(RETIRE_AGE);
for id = 1:num_id
    for t = 1:RETIRE_AGE
        sum_income_age(t) = sum_income_age(t) + income_sim(id,t);
    end
end
avg_income_age = sum_income_age./num_id; 

% Wage/Income by career
fprintf("average income by career\n");
sum_income_career = zeros(CAREER_CHOICES,1);
sum_incomesq_career = zeros(CAREER_CHOICES,1);

for id = 1:num_id
    for t = 1:7  %duration: from 20-50
        car = career_sim(id,t);
        sum_income_career(car) = sum_income_career(car) + income_sim(id,t);
    end
end
avg_income_career = sum_income_career./sum_career;

for id = 1:num_id
    for t=1:8 % from 20-55
        car = career_sim(id,t);
        sum_incomesq_career(car) = sum_incomesq_career(car) + power(income_sim(id,t) - avg_income_career(car), 2);
    end
end
sd_income_career = sqrt(sum_incomesq_career./(sum_career - 1));

% Wage/Income by edu, age and career
fprintf("average income by edu, age, and career\n");
fprintf("average income by edu, age, and career\n");
fprintf("edu, age, worker, non-incorporated, incorporated\n");
sum_income = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);

for id = 1:num_id
    for t = 1:RETIRE_AGE
        e = edu_sim(id);
        car = career_sim(id,t);
        sum_income(e,t,car) = sum_income(e,t,car) + income_sim(id,t);
    end
end
avg_income = sum_income./sum_career_ageedu;

% Wage/Income by edu and career
fprintf("average income by edu and career\n");
sum_income_educar = zeros(EDU,CAREER_CHOICES);

for id = 1:num_id
    for t = 1:RETIRE_AGE - 1
        e = edu_sim(id);
        car = career_sim(id,t);
        sum_income_educar(e,car) = sum_income_educar(e,car) + income_sim(id,t);
    end
end
avg_income_educar = sum_income_educar./sum_career_edu;

% Transition matrix
fprintf("transition matrix\n");
sum_career_lag = zeros(CAREER_CHOICES,1);
sum_career_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
sum_income_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
income_lag_vector = zeros(CAREER_CHOICES, CAREER_CHOICES); 
income_vector = zeros(CAREER_CHOICES, CAREER_CHOICES);

for id = 1:num_id
    for t = 3:8 
        car_lag = career_sim(id,t-1);
        car = career_sim(id,t);
        sum_career_tran(car_lag,car) = sum_career_tran(car_lag,car) + 1;
        sum_career_lag(car_lag) = sum_career_lag(car_lag) + 1;
        sum_income_tran(car_lag,car) = sum_income_tran(car_lag,car) + income_sim(id,t-1);
        income_lag_vector(car_lag,car) = income_sim(id,t-1);
        income_vector(car_lag,car) = income_sim(id,t);
    end
end
frac_career_tran = sum_career_tran./sum_career_lag;
avg_income_tran = sum_income_tran./sum_career_tran;

% College premium
disp("College premium");
premium_ec_sum = 0;
premium_ec_sum2 = 0;
ec_count = 0;
% period 1
for id=1:num_id
    if (edu_sim(id) == 2)
        iter = 1;
        t = 1;
        a_em = abi_em_sim(id);
        a_ub = abi_ub_sim(id);
        a_ib = abi_ib_sim(id);
        k = k_idx_sim(id,1);
        k_fam = k_fam_idx_sim(id);
        e = 1;
        exp_bus = 1;    
        
        % utility amount
        ec_count = ec_count + 1;
        premium_ec_uti = value_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam) - value_non_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam);
        premium_ec_sum = premium_ec_sum + premium_ec_uti;

        % dollar amount
        uti_k = value_worker_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k+1) - value_worker_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k);
        premium_ec_dol = CAP_STEP * premium_ec_uti ./ uti_k;
        premium_ec_sum2 = premium_ec_sum2 + premium_ec_dol;
    end
end
premium_ec_avg = premium_ec_sum ./ ec_count;
fprintf("premium_ec_avg = %f\n",premium_ec_avg);
premium_ec_avg2 = premium_ec_sum2 ./ ec_count;
fprintf("premium_ec_avg2 = %f\n",premium_ec_avg2);

total_income_elite = zeros(num_id);
total_income_nonelite = zeros(num_id);
k_idx_sim_nonelite = zeros(num_id,AGE + 1);
career_sim_nonelite = zeros(num_id,AGE);
income_sim_nonelite = zeros(num_id,AGE);
exp_bus_sim_nonelite = zeros(num_id,AGE);
total_income_diff = 0;

for id = 1:num_id
    if edu_sim(id) == 2
      iter = 1;
      a_em = abi_em_sim(id);
      a_ub = abi_ub_sim(id);
      a_ib = abi_ib_sim(id);
      k = k_idx_sim(id,1);
      k_fam = k_fam_idx_sim(id);
      total_income_elite(id) = - ELITE_TUITION + financial_aid_elite(k_fam, a_em);
      total_income_nonelite(id) = - COLLEGE_TUITION + financial_aid_nonelite(k_fam, a_em);
      schoolfuture1 = zeros(CAP);
      for k_next_idx = 1:CAP
          for consump_n_next = 1:SHOCK
              for consump_i_next = 1:SHOCK
                  schoolfuture1(k_next_idx) = schoolfuture1(k_next_idx) + BETA * value(iter,1,1,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
              end
          end
      end
      value_college = -9999;
      capitalchoice_college = -1;
      k_college = (1 + R) * (capital(k) - COLLEGE_TUITION + financial_aid_nonelite(k_fam, a_em));
      for k_next = CAP_MIN:CAP_STEP:min(k_college, capital(CAP))
          value_tmp = util(k_college - k_next) + CONSUMPTION_COLLEGE + schoolfuture1((k_next - CAP_MIN)/CAP_STEP+1);
          if (value_college < value_tmp)
              value_college = value_tmp;
              capitalchoice_college = (k_next - CAP_MIN)/CAP_STEP+1;
          end
      end
      income_sim_nonelite(id,1) = 0;
      exp_bus_sim_nonelite(id,1) = 0;
      k_idx_sim_nonelite(id,2) = capitalchoice_college;
      for t = 1:RETIRE_AGE
          e = 1;
          exp_bus = exp_bus_sim_nonelite(id,t-1);
          k = k_idx_sim_nonelite(id,t);
          consump_n = consump_n_sim(id,t);
          consump_i = consump_i_sim(id,t);
          shock_w = shock_w_sim(id,t);
          shock_n = shock_n_sim(id,t);
          shock_i = shock_i_sim(id,t);

          career_sim_nonelite(id,t) = career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i);
          if (career_sim_nonelite(id,t) == 0)
              income_sim_nonelite(id,t) = W * h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
              exp_bus_sim_nonelite(id,t) = 0;
              k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
          elseif (career_sim_nonelite(id,t) == 1)
              h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
              i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
              invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
              % liquidity constraint
              if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                  invest_non_incorp = (1 + LAMBDA_E) * capital(k);
              end
              income_sim_nonelite(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
              exp_bus_sim_nonelite(id,t) = 0;
              if (t == RETIRE_AGE)
                  exp_bus_sim_nonelite(id,t) = 1;
              end
              k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
          else
              h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
              i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
              invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
              % liquidity constraint
              if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                  invest_incorp = (1 + LAMBDA_E) * capital(k);
              end
              income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST - invest_incorp * (1 + R);
              if (exp_bus == 1)
                  income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
              end
              exp_bus_sim_nonelite(id,t) = 1;
              if (t == RETIRE_AGE)
                  exp_bus_sim_nonelite(id,t) = 2;
              end
              k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
          end
          total_income_elite(id) = total_income_elite(id) + power(BETA, t) * income_sim(id,t);
          total_income_nonelite(id) = total_income_nonelite(id) + power(BETA, t) * income_sim_nonelite(id,t);
      end
    total_income_diff = total_income_diff + total_income_elite(id) - total_income_nonelite(id);
    end
end
fprintf("elite college premium (elite college graduates) =%f", total_income_diff ./ sum_edu(3));
total_income_diff = 0;
for id = 1:num_id
    if edu_sim(id) == 1
        iter = 1;
        a_em = abi_em_sim(id);
        a_ub = abi_ub_sim(id);
        a_ib = abi_ib_sim(id);
        k = k_idx_sim(id,1);
        k_fam = k_fam_idx_sim(id);

        total_income_elite(id) = - ELITE_TUITION + financial_aid_elite(k_fam, a_em);
        total_income_nonelite(id) = - COLLEGE_TUITION + financial_aid_nonelite(k_fam, a_em);

        schoolfuture2 = zeros(CAP);
        for k_next_idx = 1:CAP
            for consump_n_next = 1:SHOCK
                for consump_i_next = 1:SHOCK
                    schoolfuture2(k_next_idx) = schoolfuture2(k_next_idx) + BETA * value(iter,1,2,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                end
            end
        end

        value_elitecollege = -9999;
        capitalchoice_elitecollege = -1;
        k_elitecollege = (1 + R) * (capital(k) - ELITE_TUITION + financial_aid_elite(k_fam, a_em));
        for k_next = CAP_MIN:CAP_STEP:min(k_elitecollege, capital(CAP))
            value_tmp = util(k_elitecollege - k_next) + CONSUMPTION_ELITE + schoolfuture2((k_next - CAP_MIN)/CAP_STEP);
            if (value_elitecollege < value_tmp)
                value_elitecollege = value_tmp;
                capitalchoice_elitecollege = (k_next - CAP_MIN)/CAP_STEP+1;
            end
        end
        if (capitalchoice_elitecollege == -1)
            capitalchoice_elitecollege = 0;
        end

        income_sim_nonelite(id,1)= 0;
        exp_bus_sim_nonelite(id,1) = 0;
        k_idx_sim_nonelite(id,2) = capitalchoice_elitecollege;
        
        for t = 1:RETIRE_AGE
            e = 2;
            exp_bus = exp_bus_sim_nonelite(id,t-1);
            k = k_idx_sim_nonelite(id,t);
            consump_n = consump_n_sim(id,t);
            consump_i = consump_i_sim(id,t);
            shock_w = shock_w_sim(id,t);
            shock_n = shock_n_sim(id,t);
            shock_i = shock_i_sim(id,t);

            if (k > 6)
                career_sim_nonelite(id,t) = career(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i);
            else
                career_sim_nonelite(id,t) = 0;
            end
            if (career_sim_nonelite(id,t) == 0)
                income_sim_nonelite(id,t) = W * h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
                exp_bus_sim_nonelite(id,t) = 0;
                k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
            elseif (career_sim_nonelite(id,t) == 1)
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                % liquidity constraint
                if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim_nonelite(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
                exp_bus_sim_nonelite(id,t) = 0;
                if (t == RETIRE_AGE)
                    exp_bus_sim_nonelite(id,t) = 1;
                end
                k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
            else
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                % liquidity constraint
                if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim_nonelite(id,t)  = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST - invest_incorp * (1 + R);
                if (exp_bus == 1)
                    income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
                end
                exp_bus_sim_nonelite(id,t) = 1;
                if (t == RETIRE_AGE)
                    exp_bus_sim_nonelite(id,t) = 2;
                end
                k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
            end
            if (k_idx_sim_nonelite(id,t+1) < 1)
              k_idx_sim_nonelite(id,t+1) = 1;
            end
            total_income_elite(id) = total_income_elite(id) + power(BETA, t) * income_sim_nonelite(id,t);
            total_income_nonelite(id) = total_income_nonelite(id) + power(BETA, t) * income_sim(id,t);
        end
      total_income_diff = total_income_diff + total_income_elite(id) - total_income_nonelite(id);
      % fout_c2 << id << ", " << total_income_elite[id] << ", " << total_income_nonelite[id] << ", " << total_income_diff << std::endl;
    end
end
% cout << "elite college premium (non-elite college graduates) =" << total_income_diff/sum_edu[1] << std::endl;
% fout_c2 << "elite college premium (non-elite college graduates) =" << total_income_diff/sum_edu[1] << std::endl;
fprintf(id, age, abi_em, abi_ub, abi_ib, own_wealth, fam_wealth, edu, career, income);

% Age of starting business
fprintf("Average age of starting business by edu");
num_age_bus = zeros(2,4);
sum_age_bus = zeros(2,4);
avg_age_bus = zeros(2,4);

for id = 1:num_id
    flag_bus_ub = 0;
    flag_bus_ib = 0;
    edu = edu_sim(id);
    for t = 1:RETIRE_AGE - 1
        if (flag_bus_ub == 0 && career_sim(id,t) == 1)
            flag_bus_ub = 1;
            num_age_bus(1,edu) = num_age_bus(1,edu) + 1;
            sum_age_bus(1,edu) = sum_age_bus(1,edu) + 20 + t*5;
            num_age_bus(1,3) = num_age_bus(1,3) + 1;
            sum_age_bus(1,3) = sum_age_bus(1,3) + 20 + t*5;
        elseif (flag_bus_ib == 0 && career_sim(id,t) == 2)
            flag_bus_ib = 1;
            num_age_bus(2,edu) = num_age_bus(2,edu) + 1;
            sum_age_bus(2,edu) = sum_age_bus(2,edu) + 20 + t*5;
            num_age_bus(2,4) = num_age_bus(2,4) + 1;
            sum_age_bus(2,4) = sum_age_bus(2,4) + 20 + t*5;
        end
    end
end
% Duation of business
disp("Duration of business by edu");
num_dur_bus = zeros(2,4);
sum_dur_bus = zeros(2,4);
avg_dur_bus = zeros(2,4);

for id=1:num_id
    flag_ub = 0;
    flag_ib = 0;
    edu = edu_sim(id);
    for t = 1:RETIRE_AGE
        if (flag_ub == 0 && career_sim(id,t) == 1)
            flag_ub = 1;
            dur_ub = 1;
        elseif (flag_ub == 1 && career_sim(id,t) == 1)
            dur_ub = dur_ub + 1;
            if (t == RETIRE_AGE - 1)
                num_dur_bus(1,edu) = num_dur_bus(1,edu) + 1;
                sum_dur_bus(1,edu) = sum_dur_bus(1,edu) + dur_ub * 5;
                num_dur_bus(1,4) = num_dur_bus(1,4) + 1;
                sum_dur_bus(1,4) = sum_dur_bus(1,4) + dur_ub * 5;
            end
        elseif (flag_ub == 1 && career_sim(id,t) ~= 1)
            flag_ub = 0;
            num_dur_bus(1,edu) = num_dur_bus(1,edu) + 1;
            sum_dur_bus(1,edu) = sum_dur_bus(1,edu) + dur_ub * 5;
            num_dur_bus(1,4) = num_dur_bus(1,4) + 1;
            sum_dur_bus(1,4) = sum_dur_bus(1,4) + dur_ub * 5;
        elseif (flag_ib == 0 && career_sim(id,t) == 2)
            flag_ib = 1;
            dur_ib = 1;
        elseif (flag_ib == 1 && career_sim(id,t) == 2)
            dur_ib = dur_ib + 1;
            if (t == RETIRE_AGE - 2)
                num_dur_bus(1,edu) = num_dur_bus(1,edu) + 1;
                sum_dur_bus(1,edu) = sum_dur_bus(1,edu) + dur_ib * 5;
                num_dur_bus(1,3) = num_dur_bus(1,3) + 1;
                sum_dur_bus(1,3) = sum_dur_bus(1,3) + dur_ib * 5;
            elseif (flag_ib == 1 && career_sim(id,t) ~= 2)
                flag_ib = 0;
                num_dur_bus(2,edu) = num_dur_bus(2,edu) + 1;
                sum_dur_bus(2,edu) = sum_dur_bus(2,edu) + dur_ib * 5;
                num_dur_bus(2,4) = num_dur_bus(2,4) + 1;
                sum_dur_bus(2,4) = sum_dur_bus(2,4) + dur_ib * 5;
            end
        end
    end
end

for bus = 1:2
    for edu=1:4
        avg_dur_bus(bus,edu) = sum_dur_bus(bus,edu)/num_dur_bus(bus,edu);
        fprintf("avg_dur_bus[%d]=%f\n",edu,avg_dur_bus(bus,edu));
    end
end
disp("finish moment calculation");