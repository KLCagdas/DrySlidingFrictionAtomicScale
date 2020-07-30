% time and time steps
N = 100;
Nstep = 1000;

% positions
zi0 = 1;
x0 = [1:N]*3;   z0 = [zi0 .* ones(1, N)];
xt = x0(1);   xts = [xt, xt-xt/2];   xD = xt - (x0(1) - x0(2));
x = x0;   z = z0;  xx0 = [x(1, 3:end), x(1, end)+1, x(1, end)+2];   xx = xx0;    # change the name xx
zt = [1 2 3 4 5]; 
global ztBar
ztBar = 0;
global ztBars
ztBars = [0];

% velocities
vx0 = zeros(1, N);    vz0 = zeros(1, N);
vtx0 = 0;   vD = 4;
vx = vx0;   vz = vz0;   vtx = vtx0;

% equations 
deriv_LenJones = 0; 
global deriv_ztBar  
deriv_ztBar = 0;
sys_pot = 0;

% attributes
E = 0.84; s = 2.56;
kT = 0.2; kX = 4.3; kpX = 5.8; kZ = 10;
dt = 0.01;
mt = 0.1; m = 0.1;
cons_LJ = 24*E*s^6;
FN = 0.2;
lenzt = length(zt);
FNs = zeros(1, lenzt);