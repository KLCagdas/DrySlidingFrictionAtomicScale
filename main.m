# time and time steps
N = 100;
Nstep = 1000;

# positions
zi0 = 1;
x0 = 0:1:N;   z0 = zi0 .* ones(0, N);
xt0 = 0:1:N;   zt0 = 0:1:N;
x = x0;   z = z0;   xt = xt0;   zt = zt0;

# velocities
vx0 = zeros(1, N);    vz0 = zeros(1, N);
vtx0 = zeros(1, N);    vtz0 = zeros(1, N);    vD0 = zeros(1, N);
vx = vx0;   vz = vz0;   vtx = vtx0;   vtz = vtz0;   vD = vD0;

# equations
LenJones = 0; eq2 = 0; eq3 = 0; vT = 0; r = 0; deriv_LenJones = 0; perp_force = 0;  FN = 0;

# attributes
E = 0.84; s = 2.56;
kT = 2; kX = 3; kpX = 4; kZ = 5;
dt = 0.1;
mt = 0.1; mi = 0.1;
cons_LJ = 24*E*s^6;

for n=1:Nstep
  dx = x - xt;    dz = z - zt;
  rhoSqr = (dx.^2 + dz.^2);
  deriv_LenJones = cons_LJ * (rhoSqr.^-4 - 2 * s^6 * rhoSqr.^-7);
    
  # force updates
  Fix = -1 * (dx .* deriv_LenJones + kX * (x - x0) + 2*kpX * x - kX * (x(1, 3:end) + x(1, 1:end)));    # index problem; add the required values manually
  Fiz = -1 * (dz .* deriv_LenJones + kZ * (z * 2 - (z0 + 1)));
  Ftx = -1 * (dz .* deriv_LenJones + kT * (xt - xD));    # change the deriv_LenJones
  
  # velocity updates
  vx += Fix * dt/mi;   Vz += Fiz * dt/mi;   vtx += Ftx * dt/mt;
  
  # position updates
  x += vx * dt;   z += vz * dt;
  xt += vtx * dt;   zt += vtz * dt;
end

function getForce()
  for i=1:N
    rhoSqr = ((x(i) - xt(i))^2 + (z(i) - zt(i))^2);
    FN += (z(i) - zt(i)) * (rhoSqr^-4 - 2 * s^6 * rhoSqr^-7);
  end
  
  FN *= cEq;
  
end

