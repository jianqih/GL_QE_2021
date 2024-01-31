function[val] = h_ib(abi,e)
    mu_ib = [0,0.56,0.28]; %incorporated
    val = exp(mu_ib(e))*abi;
end