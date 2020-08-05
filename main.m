variables   # call the file variables.m 

for n=1:Nstep
  [ztBar, deriv_ztBar] = find_ztBar(lenzt, x, xt, xts, z, zt, cons_LJ, s, ztBars, FN);
  dx = x - xt;   dz = z - ztBar;
  rhoSqr = (dx.^2 + dz.^2);   # it is not square root
  deriv_LenJones = cons_LJ .* (rhoSqr.^-4 - 2*s^6 * rhoSqr.^-7);    # for efficiency, x or z are not included, so that it can be used in all the forces
  
  % force updates
  Fx = -1 * (dx .* deriv_LenJones + kX * (x - x0) + 2*kpX * x - kX * (xx + x(1, 1:end)));
  Fz = -1 * (dz .* deriv_LenJones + kZ * (z - z0));
  Ftx = -1 * ((sum(dz .* deriv_LenJones) + FN * deriv_ztBar + kT * (xt - xD)));
  Fx(1) = Fx(length(Fx)) = 0;
  Fz(1) = Fz(length(Fz)) = 0;
  
  % plot
  plot(x, z, 'ko', "MarkerFaceColor", "b", "MarkerSize", 8)
  pause(0.1)
 
  % velocity updates
  vx = vx + Fx * dt/m;   vz = vz + Fz * dt/m;
  vxx = vxx + Fx * dt/m;    vtx = vtx + Ftx * dt/mt;        # maybe a special force Fxx should update the velocity vxx
  
  % position updates
  x = x + vx * dt;   xx = xx + vxx * dt;   z = z + vz * dt;
  xt = xt + vtx * dt;   xts = [xt xts];   xD = xD + vD * dt;
  
  % save potential of the whole system(V_T)
  sys_pot = [sys_pot, sum(4*E* sum(s^12 * rhoSqr.^-6 - s^6 * rhoSqr.^-3)) + kT/2 * (xt - xD).^2 + kX/2 * sum((x - x0).^2) + kpX/2 * ((xx - xx0) - (x - x0)).^2 + kZ/2 * (z - z0).^2 + FN * ztBar];
end
