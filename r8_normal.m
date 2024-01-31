function[randomNumber] = r8_normal(mu,sigma) % return a scalar
    randomNumber = mu + sigma*randn();
end