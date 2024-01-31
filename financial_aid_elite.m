function[fee] = financial_aid_elite(k_fam,a_em)
    if a_em >= 5
        % # higher SAT, above 1/3
        fee = 20224-32.5*k_fam+6675;
    elseif a_em >= 2
        fee = 20224-32.5*k_fam;
    else
        fee = 20224-32.5*k_fam-7432;
    end
end