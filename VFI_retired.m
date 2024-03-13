%% Guo_Leung_2021_QE 
clc
clear
close all

parameters0220;
% declare the ability grid

% period during retirement
value_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); % 9 dims for next period.
career_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); % 9 dims
capitalchoice_retire = ones(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK,SHOCK); %9 dims

% period during working
value = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 9 dims == Omega'
career = ones(RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE,CAP,SHOCK,SHOCK); % 9 dims == Omega'
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

iter = 0;
err = 1;

tol_r = 0.000001;
K = 1;

vold = util();
% three types of choice.
futurevalue_retire(k_next_idx) = futurevalue_retire(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
% futurevalue_retire1(k_next_idx) = futurevalue_retire1(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,2,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
% futurevalue_retire2(k_next_idx) = futurevalue_retire2(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,3,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);

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
                                        i_scaler_non_incorp = PROD_N * h_ub(a_ub, e) * exp(shock_grid_n(shock_n)) * power(h_em_,RHO_UB);
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
                                            for k_next = 0:CAP_STEP:round(min(k_non_incorp, capital(CAP)))
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
                                            for k_next=0:CAP_STEP:round(min(k_incorp,capital(CAP)))
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
                                                value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                                                career_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                                                for shock_w = 1:SHOCK
                                                    capitalchoice_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
                                                end
                                            elseif value_non_incorp + consump_grid_n(consump_n) >= value_incorp + consump_grid_i(consump_i)
                                                value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_non_incorp + consump_grid_n(consump_n);
                                                career_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = NON_INCORP;
                                                for shock_n = 1:SHOCK
                                                    capitalchoice_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n) = capitalchoice_non_incorp(shock_n);
                                                end
                                            else
                                                value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_incorp + consump_grid_i(consump_i);
                                                career_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = INCORP;
                                                for shock_i = 1:SHOCK
                                                    capitalchoice_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i) = capitalchoice_incorp(shock_i);
                                                end
                                            end
                                        else % worker only can choose to retire.
                                            value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                                            career_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                                            for shock_w = 1:SHOCK
                                                capitalchoice_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
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
end