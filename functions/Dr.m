function dr = Dr(t)
    % relative diffusion coefficient
    a1 = 2.414055334520468;
    b1 = 0.1607303938434139;
    c1 = 0.7335647561434824;
    dr = (a1./(1 + b1.*t).^c1);
end