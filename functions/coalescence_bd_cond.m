function [pl,ql,pr,qr] = coalescence_bd_cond(xl, ul, xr, ur, t, sigmat)
    pl = 0;
    ql = 1;
    pr = 0;
    qr = 1;
%     pl = ul;
%     ql = (xl + sigmat(t));
%     pr = ur;
%     qr = (xr + sigmat(t));
end