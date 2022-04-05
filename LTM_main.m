% AE 414 Low Thrust Maneuver Project
close all;clear;clc;
 
% Constants
g0 = 9.81;
rEarth = 6378e3;
muEarth = 3.986e14;

%% LTM for 2 Days
% Given spacecraft/orbit data
r0 = 6698e3;
v = 500e-5;
t0 = 0;
tBurn = 172800; % 2 days => s
vOrbit0 = sqrt(muEarth/r0);

% Graviational accleration function 
g = @(r) g0*(r0/r)^2;

% ODE initial conditions
IC = [1;0;1;0]; %rho0 A0 B0 theta0
nPts = 5000;
tSpan = linspace(t0,tBurn,nPts);

% ode45 function to find rho
[t,y] = ode45(@(t,y) ltmOdeSolver(t,y,r0,v),tSpan,IC);   % y = [rho; A; B; theta]

% Plot spacecraft orbit
xPlot = r0*y(:,1).*cos(y(:,4));
yPlot = r0*y(:,1).*sin(y(:,4));

figure;
plot(xPlot/10^3,yPlot/10^3);
grid on;
axis equal;
title('Spacecraft Orbit Over Two Days');
xlabel('x [km]');
ylabel('y [km]');

% Convert time (t) to normalized time (tau)
tau = sqrt(g0/r0)*t;

% Calculate velocity (normalized)
dtdT = sqrt(r0/g0);
dPdT = y(:,2)*dtdT;
dAdT = ((y(:,3).^2-y(:,1))./y(:,1).^3)*dtdT;
vOrbitNorm0 = vOrbit0*tau;
uNorm = vOrbitNorm0 + sqrt(r0*g0)*sqrt(dPdT.^2 + (y(:,1).*(dAdT+(1./y(:,1)))));

% Find minimum velocity and dimensional time (in hours)
minVel = min(uNorm);
minVelTau = find(uNorm==minVel);
minVelTime = minVelTau*sqrt(r0/g0);

% Plot velocity vs normalized time
figure;
plot(tau,uNorm);
grid on;
title('Normalized Velocity of Spacecraft');
xlabel('Normalized Time');
ylabel('Normalized Velocity');


close all;


%% Spacecraft orbit transfer
% Given spacecraft/orbit parameters
v = 2.7e-5;
rGSO = 35786e3;
tSpan = linspace(0,3e7,nPts*10);

% Calculate rho until r0*rho = rGSO
opts = odeset('Events',@(t,y) ltmOdeEventHandler(t,y,r0,rGSO));
[t,y,te,ye,ie] = ode45(@(t,y) ltmOdeSolver(t,y,r0,v),tSpan,IC,opts); % y = [rho; A; B; theta]
transferTime = te/86400; % time to reach orbit in days
fprintf('Time to reach GSO: %.2f days\n',transferTime);

% Calculate delta V 


% Plot orbit transfer
figure;
theta = linspace(0,2*pi);
xInt = r0*cos(theta);
yInt = r0*sin(theta);
xFinal = rGSO*cos(theta);
yFinal = rGSO*sin(theta);
xTransfer = r0*y(:,1).*cos(y(:,4));
yTransfer = r0*y(:,1).*sin(y(:,4));
plot(xInt/1e3,yInt/1e3,'g',xFinal/1e3,yFinal/1e3,'r');
hold on;
plot(xTransfer/1e3,yTransfer/1e3,'k');
grid on;
axis equal;
title('LTM Transfer from LEO to GSO');
xlabel('x [km]');
ylabel('y [km]');




