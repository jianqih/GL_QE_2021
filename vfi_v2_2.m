%% Guo_Leung_2021_QE 
clc
clear
close all

parameters0220;
% declare the ability grid

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

iter = 0;
err = 1;

tol_r = 0.000001;
K = 1;

while abs(err)>tol_r
    iter = iter + 1;
    R = PROD_EM*alpha*K.^(alpha-1);
    W = PROD_EM*(1-alpha)*K.^(alpha);
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
                                % capitalchoice_worker = -1;
                                k_worker = capital(k) * (1 + R) + P * W * h_em(a_em, e) * exp((4 * ALPHA_1 + 22.7 * ALPHA_2)); % actural holding
                                if (t == AGE)
                                    value_worker = util(k_worker);
                                    capitalchoice_worker = 5; % polciy function: capital = 0 
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
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                                        invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp .* NU_UB), 1/(NU_UB - 1));
                                        % liquidity constraint
                                        if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                                            invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) + (capital(k) - invest_non_incorp) * (1 + R);
                                        if (t == AGE)
                                            value_non_incorp_tmp = util(k_non_incorp);
                                            capitalchoice_non_incorp(shock_n) = 5;
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
                                        % capitalchoice_incorp(shock_i) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                                        invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                                        % liquidity constraint
                                        if invest_incorp > (1 + LAMBDA_E) * capital(k)
                                            invest_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST + (capital(k) - invest_incorp) * (1 + R);
                                        if (exp_bus == 2)
                                            k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) + (capital(k) - invest_incorp) * (1 + R);
                                        end
                                        if (t == AGE)
                                            value_incorp_tmp = util(k_incorp);
                                            capitalchoice_incorp(shock_i) = 5;
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
                                    k_worker = capital(k) * (1 + R) + W * h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w)));
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
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                                        invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                                        % liquidity constraint
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
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e)));
                                        i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                                        invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                                        % liquidity constraint
                                        if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                                            invest_incorp = (1 + LAMBDA_E) * capital(k);
                                        end
                                        k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST + (capital(k) - invest_incorp) * (1 + R);
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
                                        % no_college = no_college +1;
                                        schoolcapital2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = -1;
                                    else
                                        schoolchoice2(iter,a_em,a_ub,a_ib,k,k_fam,consump_c,consump_e) = COLLEGE;
                                        % college = college +1;
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
                income_sim(id,t) = W * h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));
                exp_bus_sim(id,t) = 1;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                laborw_sum = laborw_sum + income_sim(id,t);
                labor_sum = labor_sum + h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));                
                capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                capital_sum = capital_sum + capital(k_idx_sim(id,t));
                capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
            elseif (career_sim(id,1) == 2)
                h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
                h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
            career_sim(id,1) = 1;
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
                income_sim(id,t) = W * h_em(a_em, e) * exp((ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w)));
                exp_bus_sim(id,t) = 1;
                k_idx_sim(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                laborw_sum = laborw_sum + income_sim(id,t);
                labor_sum = labor_sum + h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e) + shock_grid_w(e,shock_w));                
                capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                capital_sum = capital_sum + capital(k_idx_sim(id,t));
                capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
            elseif (career_sim(id,t) == 2)
                h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
                h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
                    income_sim(id,t) = P * W * h_em(a_em, e) * exp(4 * ALPHA_1 + 22.7 * ALPHA_2);
                    exp_bus_sim(id,t) = 1;
                    k_idx_sim(id,t+1) = capitalchoice_retire(iter,t-RETIRE_AGE+1,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w);
                    capital_sum = capital_sum + capital(k_idx_sim(id,t));
                    capitalw_sum = capitalw_sum + capital(k_idx_sim(id,t));
                    capitalr_sum = capitalr_sum + R * capital(k_idx_sim(id,t));
                elseif (career_sim(id,t) == 2)
                    h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
                    h_em_ = h_em(a_em, e) * exp(ALPHA_1 * exp_em(t, e) + ALPHA_2 * exp_em(t, e) * exp_em(t, e));
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
                career_sim(id,t) = -1;
            end
        end
    end
    % **********************************Calculate production function******************************************************************** 
    fprintf("production function\n");
    K_EM = abs(capitalw_sum-capitalr_sum);
    K1 = K_EM/labor_sum;
    err = norm(K1-K);
    fprintf("wage:%f interest rate:%f\n",W,R);
    fprintf("K0:%f K1:%f err:%f\n",K,K1,err);
    K = 0.8*K1+0.2*K;
end