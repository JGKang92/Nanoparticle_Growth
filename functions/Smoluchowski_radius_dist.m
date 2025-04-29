function y = Smoluchowski_radius_dist(a, phi)
    W = (a + 1).*gamma(a + 3/2) / gamma(a + 2);
    y = 2.*W.*(W.*phi).^(2*a + 1).*exp(-(W.*phi).^2) ./ gamma(a + 1);
end