% First-order Upmind method to simulate the one-dimensional continuum 
% equations for the flux of electrons in the diode

% Reference: "Beyond the Child-Langmuir limit"
% R. E. Caflisch and M. S. Rosin

clear
clc

xmin = 0;
xmax = 1;
N = 100;
dt = 0.009;
t = 0;
tmax = 0.5;

dx = (xmax - xmin)/N;
x = xmin - dx : dx : xmax + dx;

%% Set initial Conditions
phi0 = exp(-200*(x-0.25).^2);
rho0 = exp(-200*(x-0.75).^2)+1;
v0 = exp(-200*(x-0.5).^2)+1;

phi = phi0; phiupd = phi0;
rho = rho0; rhoupd = rho0;
v = v0; vupd = v0;

%% loop through time
nsteps = tmax/dt;
for n = 1 : nsteps
    %% calculate boundary conditions
    phi(1) = phi(3); phi(N+3) = phi(N+1);
    rho(1) = rho(3); rho(N+3) = rho(N+1);
    v(1) = v(3); v(N+3) = v(N+1);

    %% Calculate the FOU scheme
    for i = 2 : N + 2
        vupd(i) = ((phi(i)-phi(i-1))/dx - v(i)*(v(i)-v(i-1))/dx)*dt+v(i);
        phiupd(i) = rho(i)*(dx)^2-phi(i-1)+2*phi(i);
        rhoupd(i) = -(v(i)*(rho(i)-rho(i-1))/dx+rho(i)*(v(i)-v(i-1))/dx) * dt + rho(i);
    end
    
    %% update parameter
    t = t + dt;
    v = vupd;
    phi = phiupd;
    rho = rhoupd;
    
    %% plot solution
    figure(1);
    plot(x, v, 'bo-', 'markerfacecolor', 'b');
    shg
    pause(dt); 
end







