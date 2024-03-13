function [H] = h_j(a_j,e,j,paramS) % abi_value,education,
    H = exp(log(a_j)+paramS.mu_return(j,e));
end
