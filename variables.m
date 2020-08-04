% time and time steps
N = 100;
Nstep = 1000;

% positions
x_space = 3;    z_space = 1;    # initialize the space between atoms and the initial height of them 
x0 = [1:N] * x_space;   z0 = [z_space * ones(1, N)];
xt = 1.5;   xts = [xt, xt/2];   xD = xt - x_space;
x = x0;   z = z0;  xx0 = [x(1, 3:end), x(1, end)+ 1*x_space, x(1, end)+ 2*x_space];   xx = xx0;    # 'xx' initializes the variable x_i+1 for the motion equation of x
zt = [1 2 3 4 5]; 
global ztBar
ztBar = 0;
global ztBars
ztBars = [0];

% velocities
vx0 = zeros(1, N);    vz0 = zeros(1, N);
vxx = zeros(1, N);    vtx0 = 0;   vD = 4;
vx = vx0;   vz = vz0;   vtx = vtx0;

% equations 
deriv_LenJones = 0; 
global deriv_ztBar  
deriv_ztBar = 0;
sys_pot = [];

% attributes
E = 0.84;   s = 2.56;
kT = 0.2;   kX = 1.7;   kpX = 1.7;    kZ = 10;
dt = 0.01;
mt = 0.1;   m = 0.1;
cons_LJ = 24*E*s^6;   # for efficiency
FN = 0.2;
lenzt = length(zt);   # for efficiency
FNs = zeros(1, lenzt);
