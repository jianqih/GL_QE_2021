function[order] = group_mapping(id,num_id,CAP)
    order = round(id/num_id*CAP+1);
    if (order > CAP) % censored
        order = CAP;
    elseif (order < 1)
        order = 1;
    end
end