function [ztBar, deriv_ztBar] = find_ztBar(lenzt, x, xt, xts, z, zt, cons_LJ, s, ztBars, FN)
  for i=1:lenzt
    dx = x - xt;    dz = z - zt(i);
    rhoSqr = (dx.^2 + dz.^2);
    deriv_LenJones = cons_LJ .* (rhoSqr.^-4 - 2*s^6 * rhoSqr.^-7);
    FNs(i) = sum(dz .* deriv_LenJones * (-1));
  end
  try
    ztBar = spline(FNs, zt, FN);
    ztBars = [ztBar, ztBars];   # last value of ztBar at the beginning, the first value at the end 
    deriv_ztBar = (ztBars(2) - ztBar) / (xts(2) - xt);    # written like this due to ztBars' structure
  catch
    printf("No result is found for ztBar expected from interpolation of FN\n");
    return
  end_try_catch 
end
