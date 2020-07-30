variables

for n=1:Nstep
  [ztBar, deriv_ztBar] = find_ztBar(lenzt, x, xt, xts, z, zt, cons_LJ, s, ztBars, FN);
  dx = x - xt;   dz = z - ztBar;
  rhoSqr = (dx.^2 + dz.^2);
  deriv_LenJones = cons_LJ .* (rhoSqr.^-4 - 2*s^6 * rhoSqr.^-7);
  
  % force updates
  Fx = -1 * (dx .* deriv_LenJones + kX * (x - x0) + 2*kpX * x - kX * (xx + x(1, 1:end)));   % try to write xx in a matrix form, but not as a new matrix
  Fz = -1 * (dz .* deriv_LenJones + kZ * (z - z0));
  Ftx = -1 * ((sum(dz .* deriv_LenJones) + FN * deriv_ztBar + kT * (xt - xD)));    # check if the derivative is correct in terms of zt_Bar
  Fx(1) = Fx(length(Fx)) = 0;
  Fz(1) = Fz(length(Fz)) = 0;
  plot(x, z, '*')
  pause(0.1)
  
  % velocity updates
  vx = vx + Fx * dt/m;   vz = vz + Fz * dt/m;     vtx = vtx + Ftx * dt/mt;
  
  % potential of the whole system (VT) update
  #sys_pot += 4*E*(s^12 / rhoSqr^6 - s^6 / rhoSqr^3) + kT/2 * (xt - xD)^2 + kX/2 * (x - x0)^2 + kpX/2 * ((xx - xx0) - (x - x0))^2 + kZ/2 * (z - z0)^2 + FN * ztBar; 
  # check whether (xx - x) is in coordination with the article
  
  % position updates
  x = x + vx * dt;   z = z + vz * dt;
  xt = xt + vtx * dt;   xts = [xt xts];   xD = xD + vD * dt;
  sys_pot = 0;
end
