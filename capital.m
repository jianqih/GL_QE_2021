function[val] = capital(k_index)
    CAP = 10; % par.n
    CAP_MIN = 0; % left terminal of grid
    CAP_MAX = 0.2;
    CAP_STEP = (CAP_MAX-CAP_MIN)/(CAP-1); % step
    val = (k_index - 1)*CAP_STEP+CAP_MIN;
end