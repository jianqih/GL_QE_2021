function[u]=util(c)
    SIGMA = 1.5;
    u = (c.^(1-SIGMA))/(1-SIGMA);
end