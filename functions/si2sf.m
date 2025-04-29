function sf = si2sf(si)
    % Equation (M15)
    sf = (3 .* (sqrt(3) - 1) .* si .* (4 - si) + 12 - 8 .* sqrt(3)) ./ ...
        ((3 .* si .* (12 - 6 .* si + si .^ 2) - 16) .^ (2/3) .* pi .^ (1/3));
    sf(si<0) = 1;
end