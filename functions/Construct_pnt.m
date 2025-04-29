function pnt_guess_result = Construct_pnt(expdata)
    moments = expdata.moment_raw_value;
    timedata = expdata.timedata;
    lent = length(timedata);
    maxnum = max(expdata.maxnum);
    pnt_guess_result = zeros(lent, maxnum);
    pnt_guess = @(x, a1, a2, b0, b1, b2, b3, b4, b5, b6) 1./(sqrt(2*pi*a2.^2).*x).* ...
        exp(-(log(x) - a1).^2 / (2 * a2.^2)).* (b0 + b1.*x + b2.*x.^2 + b3.*x.^3 + ...
        b4.*x.^4 + b5.*x.^5 + b6.*x.^6);
    options = optimoptions('fsolve','Display','none', ...
        'MaxFunctionEvaluations', 1e8, 'MaxIterations', 1e8', 'Algorithm','trust-region', ...
        'StepTolerance', 1e-15);
    for idx=1:lent
        momentdata = [1, moments(idx,:)];
        fun = @(x)pnt_function(x, momentdata);
        a1init = 2*log(momentdata(2)) - 0.5*log(momentdata(3));
        a2init = sqrt(log(momentdata(3)) - 2*log(momentdata(2)));
        x0 = zeros(1,9);
        x0(1) = a1init;
        x0(2) = a2init;
        x0(3) = 1.1;
        x = fsolve(fun, x0, options);
        % x = fsolve(fun, x0);
        pnt_guess_result(idx,:) = pnt_guess(1:maxnum, x(1), x(2), x(3), x(4), ...
            x(5), x(6), x(7), x(8), x(9));
        sprintf('%d of %d\n', idx, lent);
    end
    function F = pnt_function(x, momentdata)
        a1 = x(1);
        a2 = x(2);
        b0 = x(3);
        b1 = x(4);
        b2 = x(5);
        b3 = x(6);
        b4 = x(7);
        b5 = x(8);
        b6 = x(9);
        pntmoment = @(q, a1, a2, b0, b1, b2, b3, b4, b5, b6) b0.*exp(0.5*q.*(2*a1 + a2.*a2.*q)) + ...
                                                            b1.*exp(0.5*(1+q).*(2*a1 + a2.*a2.*(1+q))) + ...
                                                            b2.*exp(0.5*(2+q).*(2*a1 + a2.*a2.*(2+q))) + ...
                                                            b3.*exp(0.5*(3+q).*(2*a1 + a2.*a2.*(3+q))) + ...
                                                            b4.*exp(0.5*(4+q).*(2*a1 + a2.*a2.*(4+q))) + ...
                                                            b5.*exp(0.5*(5+q).*(2*a1 + a2.*a2.*(5+q))) + ...
                                                            b6.*exp(0.5*(6+q).*(2*a1 + a2.*a2.*(6+q)));
        qdata = 0:8;
        F = pntmoment(qdata, a1, a2, b0, b1, b2, b3, b4, b5, b6) - momentdata(1:9);
    end
end