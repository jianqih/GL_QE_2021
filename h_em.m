function[val] = h_em(abi,e)
    mu_em = [0,0.47,0.25];
    val = exp(mu_em(e))*abi;
end