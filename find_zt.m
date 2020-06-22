# time and time steps
N = 100;
Nstep = 1000;

# positions
zi0 = 1;    ztBar = 0;
x0 = 0:1:N;   z0 = zi0 .* ones(0, N);
xt0 = 0:1:N;   zt0 = 0:1:N;
x = x0;   z = z0;   xt = xt0;   zt = [1 2 3 4 5];
lenzt = length(zt);

FN = 3.2;
valFN = ones(1, lenzt);
E = 0.84; s = 2.56;
cons_LJ = 24*E*s^6;

function find_zt(xt)    # add parameters required
  for i=1:lenzt   # is it more efficient that much?
    dx = x - xt;    dz = z - zt(i);   # xt and zt must be values, not matrices
    rhoSqr = (dx.^2 + dz.^2);
    deriv_LenJones = cons_LJ * (rhoSqr.^-4 - 2 * s^6 * rhoSqr.^-7);    # kontolr et
    valFN(i) = dz .* deriv_LenJones;    # "sum function!!!"
  endfor

    # spline
  ztBar = spline(valFN, zt, FN);    # error function
endfunction
