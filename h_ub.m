function[val] = h_ub(abi,e)
    mu_ub = [0,0.35,0.20];
    val = exp(mu_ub(e))*abi;
end