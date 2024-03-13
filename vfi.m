function [Value,Polica, Polciy]= vfi()
    % Initialization
    value_retire = zeros(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); %# 10 dims
    career_retire = zeros(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK); % 10 dims
    capitalchoice_retire = zeros(AGE-RETIRE_AGE,EDU,ABI,ABI,ABI,EXP_BUS_TYPE_RETIRE,CAP,SHOCK,SHOCK,SHOCK);
    % state variable and control variable
    [e, a_em, a_ub, a_ib, exp_bus, k_next_idx, k]=ndgrid(grid.e, grid.a, grid.a, grid.a, grid.exp, grid.k, grid.k);

    for t = AGE:-1:RETIRE_AGE+1
        if t < AGE
            futurevalue_retire0 = futurevalue_retire0 + BETA * ZETA .* value_retire(e,a_em,a_ub,a_ib,1,k_next_idx) * prob_consump_n * prob_consump_i;
            futurevalue_retire1 = futurevalue_retire1 + BETA * ZETA .* value_retire(e,a_em,a_ub,a_ib,2,k_next_idx) .* prob_consump_n .* prob_consump_i;
            futurevalue_retire2 = futurevalue_retire2 + BETA * ZETA .* value_retire(e,a_em,a_ub,a_ib,3,k_next_idx) .* prob_consump_n .* prob_consump_i;
        end
        % worker's decision
        value_worker = -9999*grid.k;
        kpol = -1;
        k_worker = capital(grid.k) * (1 + R) + P * W * h_em * exp((4 * ALPHA_1 + 22.7 * ALPHA_2)); % actural holding
        if (t == AGE)
            value_worker = util(k_worker);
            kpol = 1; % polciy function: capital = 0 
        else
            k_next = interp1(k_worker,grid.k,grid.k,'linear','extrap'); % search for optimal k'
                value_tmp = util(k_worker - k_next) + futurevalue_retire0((k_next - CAP_MIN)/CAP_STEP+1);
                % disp((k_next - CAP_MIN)/CAP_STEP+1);
                if (value_worker < value_tmp)
                    value_worker = value_tmp;
                    kpol = k_next;
                end
            end
            % non-incorporated business owner decision
            value_non_incorp = 0;
            capitalchoice_non_incorp = ones(SHOCK,1);
            % only able to run a business when there is no debt
            if (capital(grid.k) >= 0)
                value_non_incorp_tmp = -9999*grid.k;
                % capitalchoice_non_incorp(shock_n) = -1;
                h_em_ = h_em .* exp((ALPHA_1 .* exp_em(t,:) + ALPHA_2 .* exp_em(t,:) .* exp_em(t,:)));
                i_scaler_non_incorp = PROD .* h_ub .* exp(shock_grid_n) * power(h_em_, RHO_UB);
                invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp .* NU_UB), 1/(NU_UB - 1));
                % liquidity constraint
                invest_non_incorp(invest_non_incorp > (1 + LAMBDA_E) * capital) = (1 + LAMBDA_E) * capital;
                k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB)- I_COST + (capital - invest_non_incorp) * (1 + R);
                if (t == AGE)
                    value_non_incorp_tmp = util(k_non_incorp);
                    capitalchoice_non_incorp(shock_n) = 1;
                else
                    k_next = interp1(k_non_incorp,grid.k,grid.k,'linear','extrap')
                    value_tmp = util(k_non_incorp - k_next) + futurevalue_retire1((k_next - CAP_MIN)/CAP_STEP+1);
                    value_non_incorp_tmp(value_non_incorp_tmp < value_tmp) = value_tmp;
                capitalchoice_non_incorp(shock_n) = (k_next - CAP_MIN)/CAP_STEP+1;
                    end
                end
                value_non_incorp = value_non_incorp + value_non_incorp_tmp * prob_shock_n(shock_n);
            end
        else
            value_non_incorp = -999;
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
                    value_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = value_worker;
                    career_retire(t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i) = WORKER;
                    for shock_w = 1:SHOCK
                        capitalchoice_retire(iter,t-RETIRE_AGE,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
                    end
                end
            end
        end
    end