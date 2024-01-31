
% //**********************************Calculate production function********************************************************************//
fprintf("production function\n");
prod_alpha = capitalr_sum/(capitalr_sum + laborw_sum);
prod_a = (R / prod_alpha) * power((labor_sum/capital_sum), prod_alpha - 1);
fprintf("prod_alpha = %f",prod_alpha,", prod_a = %f", prod_a, ", capitalr_sum = %f", capitalr_sum, ", capital_sum = %f", capital_sum, ", laborw_sum = %f", laborw_sum, ", labor_sum = %f", labor_sum);

% ************************************** Moment Calculation **************************************//
fprintf("finish simulation, start moment calculation\n");
% Proportion by edu group
fprintf("education choice\n");
sum_edu = zeros(EDU);
for id = 1:num_id
    e = edu_sim(id);
    sum_edu(e) = sum_edu(e) + 1;
end
frac_edu = sum_edu / num_id;
% Education by worker ability and wealth
% fout_est << "education choice by employee ability and wealth" << std::endl;
% fout_est << "abi_em, k , frac_hs, frac_nc, frac_ec" << std::endl;
sum_abi_wealth = zeros(ABI,6);
sum_edu_abi_wealth = zeros(ABI,6,EDU);

for id = 1:num_id
    k_group = capital(k_idx_sim(id,1))/20000;
    if k_group < 1
        k_group = 1;
    if (k_group > 5)
        k_group = 5;
    end
    e = edu_sim(id);
    sum_abi_wealth(abi_em_sim(id),k_group) = sum_abi_wealth(abi_em_sim(id),k_group) + 1;
    sum_edu_abi_wealth(abi_em_sim(id),k_group,e) = sum_edu_abi_wealth(abi_em_sim(id),k_group,e) + 1;
end

for abi_em=1:ABI
    for k = 1:6
        % fout_est << abi_em << ", " << k << ", " << (double) sum_edu_abi_wealth[abi_em][k][0] / sum_abi_wealth[abi_em][k] << ", " << (double) sum_edu_abi_wealth[abi_em][k][1] / sum_abi_wealth[abi_em][k] << ", " << (double) sum_edu_abi_wealth[abi_em][k][2] / sum_abi_wealth[abi_em][k] << ", " << sum_abi_wealth[abi_em][k] << std::endl;
    end
end

% Ability by edu group
fprintf("ability by education group\n");
sum_abi_edu_em = zeros(EDU);
frac_abi_edu_em = zeros(EDU);
sum_abi_edu_ub = zeros(EDU);
frac_abi_edu_ub = zeros(EDU);
sum_abi_edu_ib = zeros(EDU);
frac_abi_edu_ib = zeros(EDU);
for id = 1:num_id
    e = edu_sim(id);
    sum_abi_edu_em(e) = sum_abi_edu_em(e) + abi_grid_em(abi_em_sim(id));
    sum_abi_edu_ub(e) = sum_abi_edu_ub(e) + abi_grid_ub(abi_ub_sim(id));
    sum_abi_edu_ib(e) = sum_abi_edu_ib(e) + abi_grid_ib(abi_ib_sim(id));
end
for e = 1:EDU
    frac_abi_edu_em(e) = sum_abi_edu_em(e) / sum_edu(e);
    frac_abi_edu_ub(e) = sum_abi_edu_ub(e) / sum_edu(e);
    frac_abi_edu_ib(e) = sum_abi_edu_ib(e) / sum_edu(e);
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
sum_career = zeros(CAREER_CHOICES);
count_career = 0;

for id = 1:num_id
    for t=1:RETIRE_AGE - 1
        car = career_sim(id,t);
        if (car >=0)
            sum_career(car) = sum_career(car) + 1;
            count_career = count_career + 1;
        end
    end
end
frac_career = sum_career/count_career;
% Career choice by age and edu group
sum_career_ageedu = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);
frac_career_ageedu = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);

for id = 1:num_id
    for t=1:RETIRE_AGE
        e = edu_sim(id);
        car = career_sim(id,t);
        if (car >=0)
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
frac_career_edu = zeros(EDU,CAREER_CHOICES);
for id = 1:num_id
    for t = 2:RETIRE_AGE - 1
        e = edu_sim(id);
        car = career_sim(id,t);
        if (car >=0)
            sum_career_edu(e,car) = sum_career_edu(e,car) + 1;
        end
    end
end

for e = 1:EDU
    for car = 1:CAREER_CHOICES
        frac_career_edu(e,car) = sum_career_edu(e,car)/sum_edu(e)/(RETIRE_AGE - 2);
    end
end

% // Wage/Income by edu
fprintf("mean and sd of income by edu\n");
sum_income_edu = zeros(EDU);
avg_income_edu = zeros(EDU);
sum_incomesq_edu = zeros(EDU);
sd_income_edu = zeros(EDU);
for id = 1:num_id
    for t = 2:8
        e = edu_sim(id);
        sum_income_edu(e) = sum_income_edu(e) + income_sim(id,t);
    end
end
% avg_income_edu = sum_income_edu/sum_
for e = 1:EDU
    avg_income_edu(e) = sum_income_edu(e)/sum_edu(e)/(8-1);
    % fout << "avg_income_edu[" << e << "] = " << avg_income_edu[e] << std::endl;
end
for id = 1:num_id
    for t = 2:8
        e = edu_sim(id);
        sum_incomesq_edu(e) = sum_incomesq_edu(e) + power(income_sim(id,t) - avg_income_edu(e), 2);
    end
end

for e =1:EDU
    sd_income_edu(e) = sqrt(sum_incomesq_edu(e)/(sum_edu(e)*(8-1)-1));
    % fout << "sd_income_edu[" << e << "] = " << sd_income_edu[e] << std::endl;
end

% // Wage/Income by age
fprintf("average income by age\n");
sum_income_age = zeros(RETIRE_AGE);
for id = 1:num_id
    for t = 1:RETIRE_AGE
        sum_income_age(t) = sum_income_age(t) + income_sim(id,t);
    end
end
avg_income_age = sum_income_age/num_id; 

% // Wage/Income by career
fprintf("average income by career\n");
sum_income_career = zeros(CAREER_CHOICES);
avg_income_career = zeros(CAREER_CHOICES);
sum_incomesq_career = zeros(CAREER_CHOICES);
sd_income_career = zeros(CAREER_CHOICES);

for id = 1:num_id
    for t = 1:7 %duration: 7*5 = 35 
        car = career_sim(id,t);
        sum_income_career(car) = sum_income_career(car) + income_sim(id,t);
    end
end

for car=1:CAREER_CHOICES
    avg_income_career(car) = sum_income_career(car)/sum_career(car);
end

for id = 1:num_id
    for t=1:8
        car = career_sim(id,t);
        sum_incomesq_career(car) = sum_incomesq_career(car) + power(income_sim(id,t) - avg_income_career(car), 2);
    end
end

for car=1:CAREER_CHOICES
    sd_income_career(car) = sqrt(sum_incomesq_career(car)/(sum_career(car) - 1));
    % fout << "sd_income_career[" << car << "] = " << sd_income_career(car) << std::endl;
end


% // Wage/Income by edu, age and career
fprintf("average income by edu, age, and career\n");
fprintf("average income by edu, age, and career");
fprintf("edu, age, worker, non-incorporated, incorporated");
sum_income = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);
avg_income = zeros(EDU,RETIRE_AGE,CAREER_CHOICES);

for id = 1:num_id
    for t = 1:RETIRE_AGE
        e = edu_sim(id);
        car = career_sim(id,t);
        sum_income(e,t,car) = sum_income(e,t,car) + income_sim(id,t);
    end
end

for id = 1:num_id
    for t = 1:RETIRE_AGE
        for car=1:CAREER_CHOICES
            avg_income(e,t,car) = sum_income(e,t,car)/sum_career_ageedu(e,t,car);
        end
    end
end

% Wage/Income by edu and career
fprintf("average income by edu and career\n");
sum_income_educar = zeros(EDU,CAREER_CHOICES);
avg_income_educar = zeros(EDU,CAREER_CHOICES);

for id = 1:num_id
    for t = 1:RETIRE_AGE - 1
        e = edu_sim(id);
        car = career_sim(id,t);
        sum_income_educar(e,car) = sum_income_educar(e,car) + income_sim(id,t);
    end
end

for e = 1:EDU
    for car = 1:CAREER_CHOICES
        avg_income_educar(e,car) = sum_income_educar(e,car)/sum_career_edu(e,car);
    end
end

% Transition matrix
fprintf("transition matrix\n");
sum_career_lag = zeros(CAREER_CHOICES);
sum_career_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
frac_career_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
sum_income_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
avg_income_tran = zeros(CAREER_CHOICES,CAREER_CHOICES);
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


for car_lag=1:CAREER_CHOICES
    for car = 1:CAREER_CHOICES
        frac_career_tran(car_lag,car) = sum_career_tran(car_lag,car)/sum_career_lag(car_lag);
    end
end

for car_lag=1:CAREER_CHOICES
    for car = 1:CAREER_CHOICES
        avg_income_tran(car_lag,car) = sum_income_tran(car_lag,car)/sum_career_tran(car_lag,car);
    end
end

% Income correlation
fprintf("Income correlation between period t and t+1\n");
for car_lag=1:CAREER_CHOICES
    for car = 1:CAREER_CHOICES
        fout << "corr_income[" << car_lag << "][" << car << "] = " << pearsoncoeff(income_lag_vector[car_lag][car], income_vector[car_lag][car]) << std::endl;
    end
end


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

        % // utility amount
        ec_count = ec_count + 1;
        premium_ec_uti = value_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam) - value_non_elite_noshock(iter,a_em,a_ub,a_ib,k,k_fam);
        premium_ec_sum = premium_ec_sum + premium_ec_uti;

        % // dollar amount
        uti_k = value_worker_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k+1) - value_worker_noshock(iter,t,e,a_em,a_ub,a_ib,exp_bus,k);
        premium_ec_dol = CAP_STEP * premium_ec_uti / uti_k;
        premium_ec_sum2 = premium_ec_sum2 + premium_ec_dol;
    end
end
premium_ec_avg = premium_ec_sum / ec_count;
fprintf("premium_ec_avg = %f\n",premium_ec_avg);
premium_ec_avg2 = premium_ec_sum2 / ec_count;
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
      % capitalchoice_college = -1;
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
              % // liquidity constraint
              if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                  invest_non_incorp = (1 + LAMBDA_E) * capital(k);
              end
              income_sim_nonelite(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
              exp_bus_sim_nonelite(id,t) = 0;
              if (t == RETIRE_AGE - 1)
                  exp_bus_sim_nonelite(id,t) = 1;
              end
              k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
          else
              h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
              i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
              invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
              % // liquidity constraint
              if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                  invest_incorp = (1 + LAMBDA_E) * capital(k);
              end
              income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST - invest_incorp * (1 + R);
              if (exp_bus == 1)
                  income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
              end
              exp_bus_sim_nonelite(id,t) = 1;
              if (t == RETIRE_AGE - 1)
                  exp_bus_sim_nonelite(id,t) = 2;
              end
              k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
          end
          total_income_elite(id) = total_income_elite(id) + power(BETA, t) * income_sim(id,t);
          total_income_nonelite(id) = total_income_nonelite(id) + power(BETA, t) * income_sim_nonelite(id,t);
      end
    total_income_diff = total_income_diff + total_income_elite(id) - total_income_nonelite(id);
    % fout_c2 << id << ", " << total_income_elite[id] << ", " << total_income_nonelite[id] << ", " << total_income_diff << std::endl;
    end
end
fprintf("elite college premium (elite college graduates) =%f", total_income_diff/sum_edu(3));

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
        % capitalchoice_elitecollege = -1;
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
                % // liquidity constraint
                if (invest_non_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_non_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim_nonelite(id,t) = (1-DELTA) * invest_non_incorp + i_scaler_non_incorp * power(invest_non_incorp, NU_UB) - invest_non_incorp * (1 + R);
                exp_bus_sim_nonelite(id,t) = 0;
                if (t == RETIRE_AGE - 1)
                    exp_bus_sim_nonelite(id,t) = 1;
                end
                k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_n);
            else
                h_em_ = h_em(a_em, e) * exp(ALPHA_1(e) * exp_em(t, e) + ALPHA_2(e) * exp_em(t, e) * exp_em(t, e));
                i_scaler_incorp = PROD_I * h_ib(a_ib, e) * exp(shock_grid_i(shock_i)) * power(h_em_, RHO_IB);
                invest_incorp = power((DELTA+R)/(i_scaler_incorp * NU_IB), 1/(NU_IB - 1));
                % // liquidity constraint
                if (invest_incorp > (1 + LAMBDA_E) * capital(k))
                    invest_incorp = (1 + LAMBDA_E) * capital(k);
                end
                income_sim_nonelite(id,t)  = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - I_COST - invest_incorp * (1 + R);
                if (exp_bus == 1)
                    income_sim_nonelite(id,t) = (1-DELTA) * invest_incorp + i_scaler_incorp * power(invest_incorp, NU_IB) - invest_incorp * (1 + R);
                end
                exp_bus_sim_nonelite(id,t) = 1;
                if (t == RETIRE_AGE - 1)
                    exp_bus_sim_nonelite(id,t) = 2;
                end
                k_idx_sim_nonelite(id,t+1) = capitalchoice(iter,t,e,a_em,a_ub,a_ib,exp_bus,k,consump_n,consump_i,shock_i);
            end
            if (k_idx_sim_nonelite(id,t+1) < 0)
              k_idx_sim_nonelite(id,t+1) = 0;
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
for id = 1:num_id
    for t =1:RETIRE_AGE
        % fout_c1 << id << ", " << t << ", " << abi_em_sim[id] << ", " << abi_ub_sim[id] << ", " << abi_ib_sim[id] << ", " << k_idx_sim[id][0] << ", " << k_fam_idx_sim[id] << ", " <<  edu_sim[id] << ", " << career_sim[id][t] << ", " << income_sim[id][t] << std::endl;
    end
end

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

for bus=1:2
    for edu=1:4
        avg_age_bus(bus,edu) = sum_age_bus(bus,edu)/num_age_bus(bus,edu);
        % fout << "avg_age_bus[" << bus << "][" << edu << "] = " << avg_age_bus[bus][edu] << std::endl;
    end
end

% // Duation of business
disp("Duration of business by edu");
num_dur_bus = zeros(2,4);
sum_dur_bus = zeros(2,4);
avg_dur_bus = zeros(2,4);

for id=1:num_id
    flag_ub = 0;
    flag_ib = 0;
    edu = edu_sim(id);
    for t = 1:RETIRE_AGE - 1
        if (flag_ub == 0 && career_sim(id,t) == 1)
            flag_ub = 1;
            dur_ub = 1;
        elseif (flag_ub == 1 && career_sim(id,t) == 1)
            dur_ub = dur_ub + 1;
            if (t == RETIRE_AGE - 2)
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
