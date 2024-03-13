%% Guo_Leung_2021_QE 
clc
clear
close all

parameters;
% declare the ability grid

value_retire = ones(cS.AGE-cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE_RETIRE,cS.nk,cS.shock,cS.shock); %# 10 dims
career_retire = ones(cS.AGE-cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE_RETIRE,cS.nk,cS.shock,cS.shock); %# 10 dims
capitalchoice_retire = ones(cS.AGE-cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE_RETIRE,cS.nk,cS.shock,cS.shock); %# 10 dims

value = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk,cS.shock,cS.shock); % 10 dims
career = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk,cS.shock,cS.shock); % 10 dims
capitalchoice = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk,cS.shock,cS.shock,cS.shock); % 10 dims
transfer = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk,cS.shock,cS.shock); % 10 dims
value_worker_noshock = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk); % 10 dims
value_non_incorp_noshock = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk); % 10 dims
value_incorp_noshock = ones(cS.RETIRE_AGE,EDU,cS.abi,cS.abi,cS.abi,EXP_BUS_TYPE,cS.nk); % 10 dims

schoolvalue = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk,cS.shock,cS.shock); % 8 dims
schoolchoice = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk,cS.shock,cS.shock); % 8 dims
schoolcapital = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk,cS.shock,cS.shock); % 8 dims
schoolchoice2 = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk,cS.shock,cS.shock); % 8 dims
schoolcapital2 = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk,cS.shock,cS.shock); % 8 dims
value_non_elite_noshock = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk); % 8 dims
value_elite_noshock = ones(cS.abi,cS.abi,cS.abi,cS.nk,cS.nk); % 8 dims

borrow_incorp = zeros(cS.AGE,cS.nk);
borrow_non_incorp = zeros(cS.AGE,cS.nk);
K_worker = zeros(cS.AGE,cS.nk);
iter = 0;
err = 1;
minrate = cS.delta;
tol_r = 0.000001;
maxrate = (1-cS.beta)/(cS.beta);

K = 1;

while abs(err)>tol_r
    iter = iter + 1;
    R = alpha * PROD_EM * (K.^(alpha-1));
    W = (1-alpha) * PROD_EM * K.^(alpha);
    fprintf("Iteration: %d\n",iter);
    for t = cS.AGE:-1:cS.RETIRE_AGE+1
        for exp_bus=1:EXP_BUS_TYPE_RETIRE
        fprintf("Processing Retirement age\n");
        fprintf("Processing age = %d\n",age(t));
        futurevalue_retire0 = zeros(cS.nk,1);
        futurevalue_retire1 = zeros(cS.nk,1);
        futurevalue_retire2 = zeros(cS.nk,1);
        for k_next_idx = 1:cS.nk
            if (t < cS.AGE)
                for consump_n_next=1:cS.shock
                    for consump_i_next=1:cS.shock
                        futurevalue_retire0(k_next_idx) = futurevalue_retire0(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,:,:,:,:,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        futurevalue_retire1(k_next_idx) = futurevalue_retire1(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,:,:,:,:,2,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        futurevalue_retire2(k_next_idx) = futurevalue_retire2(k_next_idx) + BETA * ZETA(t-RETIRE_AGE) * value_retire(t-RETIRE_AGE,:,:,:,:,3,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                    end
                end
            end
        end
        for k = 1:cS.nk % index k 
            % worker's decision
            value_worker = -9999;
            capitalchoice_worker = -1;
            k_worker = capital(k) * (1 + R) + cS.p * W * h_em * exp((4 * ALPHA_1 + 22.7 * ALPHA_2)); % actural holding
            if (t == cS.AGE)
                value_worker = util(k_worker,cS);
                capitalchoice_worker = 5; % polciy function: capital = 0 
            else
                for k_next = 0:cS.kstep:round(min(k_worker, capital(cS.nk))) % real capital k next period
                    value_tmp = util(k_worker - k_next) + futurevalue_retire0((k_next - cS.kMin)/cS.kstep+1);
                    % disp((k_next - CAP_MIN)/CAP_STEP+1);
                    if (value_worker < value_tmp)
                        value_worker = value_tmp;
                        capitalchoice_worker = (k_next - CAP_MIN)/CAP_STEP+1;
                    end
                    K_worker(t,k) = capital(capitalchoice_worker);
                end
            end

            % non-incorporated business owner decision
            value_non_incorp = 0;
            capitalchoice_non_incorp = ones(cS.shock,1);
            % only able to run a business when there is no debt
            if (capital(k) >= 0)
                for shock_n = 1:cS.shock
                    value_non_incorp_tmp = -9999;
                    capitalchoice_non_incorp(shock_n) = -1;
                    h_em_ = h_em(a_em, :) * exp((ALPHA_1 * exp_em(t, :) + ALPHA_2 * exp_em(t, :) * exp_em(t, :)));
                    i_scaler_non_incorp = PROD_N * h_ub(a_ub, :) * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                    invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp .* NU_UB), 1/(NU_UB - 1));
                    % liquidity constraint
                    if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                        invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    borrow_non_incorp(t,k) = invest_non_incorp - capital(k);
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
            capitalchoice_incorp = ones(cS.shock,1);
            % only able to run a business when there is no debt
            if (capital(k) >= 0)
                for shock_i=1:cS.shock
                    value_incorp_tmp = -9999;
                    capitalchoice_incorp(shock_i) = -1;
                    h_em_ = h_em * exp((ALPHA_1 * exp_em(t, :) + ALPHA_2 * exp_em(t, :) * exp_em(t, :)));
                    i_scaler_incorp = PROD_I * h_ib * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                    invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1)); % From FOC
                    % liquidity constraint
                    if invest_incorp > (1 + LAMBDA_E) * capital(k)
                        invest_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    borrow_incorp(t,k) = invest_incorp - capital(k);
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

            for consump_n = 1:cS.shock
                for consump_i = 1:cS.shock
                    if (exp_bus > 1)
                        if value_worker > value_non_incorp + consump_grid_n(consump_n) && value_worker > value_incorp + consump_grid_i(consump_i)
                            value_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = value_worker;
                            career_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = WORKER;
                            for shock_w = 1:cS.shock
                                capitalchoice_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
                            end
                        elseif value_non_incorp + consump_grid_n(consump_n) >= value_incorp + consump_grid_i(consump_i)
                            value_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = value_non_incorp + consump_grid_n(consump_n);
                            career_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = NON_INCORP;
                            for shock_n = 1:SHOCK
                                capitalchoice_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i,shock_n) = capitalchoice_non_incorp(shock_n);
                            end
                        else
                            value_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = value_incorp + consump_grid_i(consump_i);
                            career_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = INCORP;
                            for shock_i = 1:cS.shock
                                capitalchoice_retire(t-RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i,shock_i) = capitalchoice_incorp(shock_i);
                            end
                        end
                    else
                        value_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = value_worker;
                        career_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i) = WORKER;
                        for shock_w = 1:SHOCK
                            capitalchoice_retire(t-cS.RETIRE_AGE,:,:,:,:,exp_bus,k,consump_n,consump_i,shock_w) = capitalchoice_worker;
                        end
                    end
                end
            end
        end
        end
    end

    for t = cS.RETIRE_AGE:-1:1
        fprintf("Processing Working age\n")
        fprintf("Processing age = %d\n",age(t))
        % pragma omp parallel for collapse(4)
        futurevalue0 = zeros(cS.nk,1);
        futurevalue1 = zeros(cS.nk,1);
        futurevalue_retire0 = zeros(cS.nk,1);
        futurevalue_retire1 = zeros(cS.nk,1);
        futurevalue_retire2 = zeros(cS.nk,1);
        for k_next_idx = 1:cS.nk
            if (t == cS.RETIRE_AGE)
                for consump_n_next = 1:cS.shock
                    for consump_i_next= 1:cS.shock
                        futurevalue_retire0(k_next_idx) = futurevalue_retire0(k_next_idx) + cS.beta * value_retire(1,:,:,:,:,:,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        futurevalue_retire1(k_next_idx) = futurevalue_retire0(k_next_idx) + cS.beta * value_retire(1,:,:,:,:,:,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        futurevalue_retire2(k_next_idx) = futurevalue_retire2(k_next_idx) + cS.beta * value_retire(1,:,:,:,:,:,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                    end
                end
            else
                for consump_n_next = 1:cS.shock
                    for consump_i_next = 1:cS.shock
                        futurevalue0(k_next_idx) = futurevalue0(k_next_idx) + cS.beta * value(t+1,:,:,:,:,:,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                        futurevalue1(k_next_idx) = futurevalue0(k_next_idx) + cS.beta * value(t+1,:,:,:,:,:,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                    end
                end
            end
        end

        for k=1:cS.nk
            % worker's decision
            value_worker = 0;
            capitalchoice_worker = ones(cS.shock,1);
            for shock_w = 1:cS.shock
                value_worker_tmp = -9999;
                capitalchoice_worker(shock_w) = -1;
                k_worker = capital(k) * (1 + R) + W * h_em * exp((ALPHA_1 * exp_em(t, :) + ALPHA_2 * exp_em(t, :) * exp_em(t, :) + shock_grid_w(:,shock_w)));
                for k_next = 0:CAP_STEP:round(min(k_worker,capital(CAP))) % capital cannot negative
                    value_tmp = util(k_worker - k_next,cS);
                    if t == cS.RETIRE_AGE
                        value_tmp = value_tmp + futurevalue_retire0((k_next - cS.kMin)/cS.kstep+1);
                    else 
                        value_tmp = value_tmp + futurevalue0((k_next - cS.kMin)/cS.kstep+1);
                    end
                    if (value_worker_tmp < value_tmp)
                        value_worker_tmp = value_tmp;
                        capitalchoice_worker(shock_w) = (k_next - cS.kMin)/cS.kstep+1;
                        K_worker(t,k) = capital(capitalchoice_worker);
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
                    h_em_ = h_em * exp((ALPHA_1 * exp_em(t, :) + ALPHA_2 * exp_em(t, :) * exp_em(t, :)));
                    i_scaler_non_incorp = cS.PROD_N * h_ub * exp(shock_grid_n(shock_n)) * power(h_em_, RHO_UB);
                    invest_non_incorp = power((DELTA+R)/(i_scaler_non_incorp * NU_UB), 1/(NU_UB - 1));
                    % // liquidity constraint
                    if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                        invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    borrow_non_incorp(t,k) = invest_non_incorp - capital(k);
                    k_non_incorp = (1 - DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) + (capital(k) - invest_non_incorp) * (1 + R);
                    for k_next = 0:CAP_STEP:min(k_non_incorp, capital(CAP))
                        value_tmp = util(k_non_incorp - k_next);
                        if (t == RETIRE_AGE)
                            value_tmp = value_tmp + futurevalue_retire1((k_next - cS.kMin)/cS.kstep+1);
                        else
                            value_tmp = value_tmp + futurevalue0((k_next - cS.kMin)/cS.kstep+1);
                        end
                        if (value_non_incorp_tmp < value_tmp)
                            value_non_incorp_tmp = value_tmp;
                            capitalchoice_non_incorp(shock_n) = (k_next - cS.kMin)/cS.kstep+1;
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
                    h_em_ = h_em * exp((ALPHA_1 * exp_em(t, :) + ALPHA_2 * exp_em(t, :) * exp_em(t, :)));
                    i_scaler_incorp = cS.PROD_I * h_ib * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                    invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                    % liquidity constraint
                    if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                        invest_incorp = (1 + LAMBDA_E) * capital(k);
                    end
                    borrow_incorp(t,k) = invest_incorp - capital(k);
                    k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST + (capital(k) - invest_incorp) * (1 + R);
                    if (exp_bus == 1)
                        k_incorp = (1 - DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) + (capital(k) - invest_incorp) * (1 + R);
                    end
                    for k_next = 0:CAP_STEP:min(k_incorp,capital(cS.nk))
                        value_tmp = util(k_incorp - k_next);
                        if (t == RETIRE_AGE)
                            value_tmp = value_tmp + futurevalue_retire2((k_next - cS.kMin)/cS.kstep+1);
                        else
                            value_tmp = value_tmp + futurevalue1((k_next - cS.kMin)/cS.kstep+1);
                        end
                        if (value_incorp_tmp < value_tmp)
                            value_incorp_tmp = value_tmp;
                            capitalchoice_incorp(shock_i) = (k_next - cS.kMin)/cS.kstep+1;
                        end
                    end
                    value_incorp = value_incorp + value_incorp_tmp * prob_shock_i(shock_i);
                end
            else
                value_incorp = -999;
            end

            
            value_worker_noshock(t,:,:,:,:,:,k) = value_worker;
            value_non_incorp_noshock(t,:,:,:,:,:,k) = value_non_incorp;
            value_incorp_noshock(t,:,:,:,:,:,k) = value_incorp;

            for consump_n = 1:cS.shock
                for consump_i = 1:cS.shock
                    if (value_worker > value_non_incorp + consump_grid_n(consump_n) && value_worker > value_incorp + consump_grid_i(consump_i))
                        value(t,:,:,:,:,:,k,consump_n,consump_i) = value_worker;
                        career(t,:,:,:,:,:,k,consump_n,consump_i) = WORKER;
                        for shock_w = 1:cS.shock
                            capitalchoice(t,:,:,:,:,:,k,consump_n,consump_i,shock_w) = capitalchoice_worker(shock_w);
                        end
                    elseif (value_non_incorp + consump_grid_n(consump_n) >= value_incorp + consump_grid_i(consump_i))
                        value(t,:,:,:,:,:,k,consump_n,consump_i) = value_non_incorp + consump_grid_n(consump_n);
                        career(t,:,:,:,:,:,k,consump_n,consump_i) = NON_INCORP;
                        for shock_n = 1:cS.shock
                            capitalchoice(t,:,:,:,:,:,k,consump_n,consump_i,shock_n) = capitalchoice_non_incorp(shock_n);
                        end
                    else
                        value(t,:,:,:,:,:,k,consump_n,consump_i) = value_incorp + consump_grid_i(consump_i);
                        career(t,:,:,:,:,:,k,consump_n,consump_i) = INCORP;
                        for shock_i = 1:cS.shock
                            capitalchoice(t,:,:,:,:,:,k,consump_n,consump_i,shock_i) = capitalchoice_incorp(shock_i);
                        end
                    end
                end
            end
        end
    end

    fprintf("start schooling decision\n");
    % schooling decision
    schoolfuture1 = zeros(cS.nk,1);
    schoolfuture2 = zeros(cS.nk,1);
    for k_next_idx=1:CAP
        for consump_n_next=1:cS.shock
            for consump_i_next=1:cS.shock
                schoolfuture1(k_next_idx) = schoolfuture1(k_next_idx) + BETA * value(2,2,:,:,:,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
                schoolfuture2(k_next_idx) = schoolfuture1(k_next_idx) + BETA * value(2,3,:,:,:,1,k_next_idx,consump_n_next,consump_i_next) * prob_consump_n(consump_n_next) * prob_consump_i(consump_i_next);
            end
        end
    end

    for k = 1:CAP
        for k_fam = 1:CAP
            % high school edu = 1
            value_highschool = 0;
            for consump_n = 1:SHOCK
                for consump_i = 1:SHOCK
                    value_highschool = value_highschool + value(1,1,:,:,:,1,k,consump_n,consump_i) * prob_consump_n(consump_n) * prob_consump_i(consump_i);
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

            value_non_elite_noshock(:,:,:,k,k_fam) = value_college;
            value_elite_noshock(:,:,:,k,k_fam) = value_elitecollege;

            for consump_c = 1:SHOCK
                for consump_e = 1:SHOCK
                    if (value_highschool > value_college + consump_grid_c(consump_c) && value_highschool > value_elitecollege + consump_grid_e(consump_e))
                        schoolchoice(:,:,:,k,k_fam,consump_c,consump_e) = NO_COLLEGE;
                        schoolvalue(:,:,:,k,k_fam,consump_c,consump_e) = value_highschool;
                        schoolcapital(:,:,:,k,k_fam,consump_c,consump_e) = -1;
                    elseif (value_college + consump_grid_c(consump_c) > value_elitecollege + consump_grid_e(consump_e))
                        schoolchoice(:,:,:,k,k_fam,consump_c,consump_e) = COLLEGE;
                        schoolvalue(:,:,:,k,k_fam,consump_c,consump_e) = value_college + consump_grid_c(consump_c);
                        schoolcapital(:,:,:,k,k_fam,consump_c,consump_e) = capitalchoice_college;
                    else
                        schoolchoice(:,:,:,k,k_fam,consump_c,consump_e) = ELITE;
			            % check if the student is admitted by elite college
                        if a_em >= 5
				            sat_g = 3;
                        elseif a_em >= 2
			 	            sat_g = 2;
                        else 
                            sat_g = 1;
                        end
                        schoolvalue(:,:,:,k,k_fam,consump_c,consump_e) = ELITE_ADMIT(sat_g) * value_elitecollege + (1-ELITE_ADMIT(sat_g)) * value_college + consump_grid_e(consump_e);
                        schoolcapital(:,:,:,k,k_fam,consump_c,consump_e) = capitalchoice_elitecollege;
                        if (value_highschool > value_college + consump_grid_c(consump_c))
                            schoolchoice2(:,:,:,k,k_fam,consump_c,consump_e) = NO_COLLEGE;
                            % no_college = no_college +1;
                            schoolcapital2(:,:,:,k,k_fam,consump_c,consump_e) = -1;
                        else
                            schoolchoice2(:,:,:,k,k_fam,consump_c,consump_e) = COLLEGE;
                            % college = college +1;
                            schoolcapital2(:,:,:,k,k_fam,consump_c,consump_e) = capitalchoice_college;
                        end
                    end
                end
            end
        end
    end
    
    K1 = sum(sum(K_worker - borrow_incorp - borrow_non_incorp));
    R1 = PROD_EM*(K1/labor_em)^(1-alpha)*alpha;

    err = R1-R;
    if err < 0
        maxrate = R;
    else
        minrate = R;
    end
    disp(['k0 = ',num2str(K0),', k1 = ',num2str(K1),', r0 = ',num2str(R),', r1 = ',num2str(R1), ', err = ', num2str(err)])
end