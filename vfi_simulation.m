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
minrate = DELTA;
tol_r = 0.000001;
maxrate = (1-BETA)/(BETA);
R = (minrate+maxrate)*0.5; %

while abs(err)>tol_r
    iter = iter + 1;
    R = 0.5*(maxrate+minrate);
    K0 = labor_em*(R/(PROD_EM*alpha))^(1/(alpha-1));
    W = PROD_EM*(K0/labor_em)^(alpha)*(1-alpha);
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
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
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
                                        capitalchoice_incorp(shock_i) = -1;
                                        h_em_ = h_em(a_em, e) * exp((ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e)));
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
    
    
    K1 = sum(capital(capitalchoice_worker) - capital(capitalchoice_incorp) - capital(capitalchoice_non_incorp));
    R1 = PROD_EM*(max(K1,0.001)/labor_em)^(1-alpha)*alpha;

    err = R1-R;
    if err < 0
        maxrate = R;
    else
        minrate = R;
    end
    disp(['k0 = ',num2str(K0),', k1 = ',num2str(K1),', r0 = ',num2str(R),', r1 = ',num2str(R1)])
end