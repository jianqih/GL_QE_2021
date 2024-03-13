function [H] = h_em(a_em,e,t,paramS)
    H = exp(log(a_em)+paramS.mu_em(e)+paramS.gamma_1*exp_em(t,e)+paramS.gamma_2*exp_em(t,e)*exp_em(t,e));
end
