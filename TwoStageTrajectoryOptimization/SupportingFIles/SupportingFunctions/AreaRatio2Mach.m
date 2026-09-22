function machNumber = AreaRatio2Mach(areaRatio, gamma)
g = gamma;
AR = areaRatio;
gp1O2gm2 = (g+1) / (2*g-2);

AM_func = @(M) (1./M) .* (2/(g+1) .* (1 + (g-1)/2 .* M.^2)).^((g+1)/(2*(g-1))) - AR;
M_guess = [1.0, 20];
machNumber = fzero(AM_func, M_guess);
end