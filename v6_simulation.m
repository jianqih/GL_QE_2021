    %% Simulation %%
    fprintf("start simulation\n");
    [paramS.tranprob_abi_em,paramS.abi_loggrid_em,paramS.prob_abi_em] = markovappr(THETA_EM,ABI_EM_STD,RANGE_ABI,ABI);
    [paramS.tranprob_abi_ub,paramS.abi_loggrid_ub,paramS.prob_abi_ub] = markovappr(THETA_UB,ABI_UB_STD,RANGE_ABI,ABI);
    [paramS.tranprob_abi_ib,paramS.abi_loggrid_ib,paramS.prob_abi_ib] = markovappr(THETA_IB,ABI_IB_STD,RANGE_ABI,ABI);    
    % simulation array
    nSim=5000;
    kHistM = zeros(nSim, J);
    cHistM = zeros(nSim, J);
    abiEMIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_em,paramS.prob_abi_em); % (index,age)
    abiUBIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_ub,paramS.prob_abi_ub);
    abiIBIdx = AbilitySimulation_olgm(nSim,J,paramS.tranprob_abi_ib,paramS.prob_abi_ib);
    k_idx_sim = round((ABI-1).* rand(nSim,J)+1); % individual capital simulation
    k_fam_idx_sim = round((ABI-1).* rand(nSim,J)+1);


    edu_sim =  ones(num_id,1);
    
    careerWsim = ones(num_id,1);
    income_sim = ones(num_id,1);
   
    % simulate first period: education decision
    fprintf("simulation first period\n");
    for t=1:T    
        for id = 1:nSim         
            a_em = abiEMIdx(id,t);
            a_ub = abiUBIdx(id,t);
            a_ib = abiIBIdx(id,t);
            k = k_idx_sim(id,t);
            k_fam = k_fam_idx_sim(id,t);        
            edu_random_number = (randn() / (RAND_MAX)); 
        
            edu_sim(id,t) = careerEdu(a_em,a_ub,a_ib,k,k_fam);
            % check if the student is admitted by elite college
            flag = 0;
            if (edu_sim(id,t) == 3)
                sat_group=1;
	            if a_em >= 5
	                sat_group = 3;
                elseif a_em >= 2
	                sat_group = 2;
                end
	            if edu_random_number > ELITE_ADMIT(sat_group)
		            edu_sim(id,t) = careerEdu(a_em,a_ub,a_ib,k,k_fam);
		            flag = 1; % admitted successfully.
                else
                    edu_sim(id,t)=2;
                end            
            end
        end
    end


    L_EM = 0;
    for id =1:nSim
        for t = 2:T-1
            for a_em =1:na
                for e =1:ne
                    L_EM=L_EM+h_em(a_em, e)*exp(ALPHA_1*exp_em(t,e)+ALPHA_2*exp_em(t,e)*exp_em(t,e));
                end
            end
        end
    end
    for id =1:nSim
        for t = 2 : tW
           for a_em = 1 : na
               for a_ub =1:na
                   for a_ib = 1:na
                        e = edu_sim(id,t);
                    % Find next period capital for each individual by interpolation
                        kHistM(idxV, t+1, 1) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,1), ...
                                                    kHistM(idxV, t,1), 'linear');
                        kHistM(idxV, t+1, 2) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,2), ...
                                                    kHistM(idxV, t,1), 'linear');
                        kHistM(idxV, t+1, 3) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,3), ...
                                                    kHistM(idxV, t,1), 'linear');                           
                   end
               end               
           end
        end
        for t=1:tR
            for a_em = 1 : na
               for a_ub =1:na
                   for a_ib = 1:na
                        e = edu_sim(id,t);
                        % Find next period capital for each individual by interpolation
                        kHistM(idxV, t+1+tW, 1) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,1), ...
                                                    kHistM(idxV, t,1), 'linear');
                        kHistM(idxV, t+1+tW, 2) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,2), ...
                                                    kHistM(idxV, t,1), 'linear');
                        kHistM(idxV, t+1+tW, 3) = interp1(kap(:), kapWopt(a_em,a_ub,a_ib,e,:,t,3), ...
                                                    kHistM(idxV, t,1), 'linear');                           
                   end
               end               
            end
        end
   end % for ie