D = .1;
L = .001; % in meters
C0 = 100; % concentration at base of column (Pa)
Cout = 0; % fixed percentage concentration leaving top of column (%)

rho = 2335;  % density of rock ( shale = 2335 )  (kg / m^3)
g = 9.8; % m^2 / s
eta = .05; % porosity (shale = 0.0 - 0.1) (unitless)
K = 10^-14; % hydraulic conductivity ( shale =  10-15 to 10-13 )  ( m / s)

% Pmax, pressure at the bottom of the column
% according to Wikipedia, fracking equipment can reach up to 100 megapascals (1000000 Pa) (15,000 psi)

Pmin = L * rho * g + 0; % minimum P at surface to avoid downward motion
Pmax = 100000000; % some maximum P at bottom of  column

timescale = 60*60*24; % size of timestep , using 1 day

% simulating nonlinear change in P
%dPdX = @(x) -2 * (Pmax - Pmin) / L^2 * x;

% linear change in P
dPdX = @(x) (Pmin - Pmax) / L;  % just the slope

% calculate darcy velocity corrected for porosity
Vd = @(x) - eta * K * timescale * ( 1 + dPdX(x) / (rho * g) );  % m / timescale

%pdefun = @(x,t,u,DuDx) deal(1/D, DuDx, - ( Vd(x) / D) * DuDx);
pdefun = @(x,t,u,DuDx) deal(1, 0, - Vd(x) * DuDx);  % darcy flow only, no diffusion/dispersion

icfun = @(x) 0;
bcfun = @(xl,ul,xr,ur,t) deal(ul - C0, 0, Cout * ur, 1);

x = linspace(0,L,10);

years = 10;
days = 365 * years;
t = linspace(0,days,years);  % time is in days now

%years = 10;
%t = linspace(0,years*365*24*60*60,years);  % time is in seconds, convert days to seconds with a 1 day step

m = 0;
sol = pdepe(m,pdefun,icfun,bcfun,x,t);
c = sol(:,:,1);
figure; surf(x,t,c);