function irf_levels = units_to_levels(irf, tc)

% Compute impulse responses for levels

if tc == 2 || tc == 5;
    irf = cumsum(irf);
elseif tc == 3 || tc == 6;
    %irf = cumsum(cumsum(irf));
    irf = cumsum(irf);
end;

irf_levels = irf;

end