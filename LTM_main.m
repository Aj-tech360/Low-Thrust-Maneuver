% Low Thrust Maneuver Project
% Ronak Amin, Benjamin Sites, Christopher Rappole 
% AE 414 – 01 
% Prof. Laksh Narayanaswami 
% April 22, 2022 

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
r = r0*y(:,1);
xPlot = r.*cos(y(:,4));
yPlot = r.*sin(y(:,4));
figure;
plot(xPlot/rEarth,yPlot/rEarth);
grid on;
axis equal;
title('Spacecraft Orbit Over Two Days');
xlabel('x [Earth Radii]');
ylabel('y [Earth Radii]');

% Convert time (t) to normalized time (tau)
tau = sqrt(g0/r0)*t;

% Calculate velocity (dimensional)
rDot = r0*y(:,2);
thetaDotSquared = muEarth./r.^3; %(1./y(:,1)).*(y(:,2).*sqrt(r0/g0) + (1./y(:,1).^2)).*sqrt(g0/r0);
uDim = vOrbit0+sqrt(rDot.^2 + thetaDotSquared.*(r.^2));

% Find minimum velocity and dimensional time (in hours)
minVel = min(uDim);
minVelTau = find(uDim==minVel);
minVelTime = (minVelTau*sqrt(r0/g0))/3600;
fprintf('The minimum velocity is %.2f km/s at %.2f hr\n',minVel/1e3,minVelTime);

% Plot velocity vs normalized time
figure;
plot(tau,uDim/1e3);
grid on;
title('Normalized Velocity of Spacecraft');
xlabel('Normalized Time');
ylabel('Dimensional Velocity [km/s]');

%{
%% Spacecraft LTM Orbit Transfer
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
theta = linspace(0,2*pi,1000);
xInt = r0*cos(theta);
yInt = r0*sin(theta);
xFinal = rGSO*cos(theta);
yFinal = rGSO*sin(theta);
xTransfer = r0*y(:,1).*cos(y(:,4));
yTransfer = r0*y(:,1).*sin(y(:,4));
plot(xInt/rEarth,yInt/rEarth,'g',xFinal/rEarth,yFinal/rEarth,'r');
hold on;
plot(xTransfer/rEarth,yTransfer/rEarth,'color','#0072BD');
grid on;
axis equal;
title('LTM Transfer from LEO to GSO');
xlabel('x [Earth Radii]');
ylabel('y [Earth Radii]');


%% Hohmann Tranfer Calculations
% Tranfer orbit calculations
rTransfer = (r0+rGSO)/2;
eTransfer = -muEarth/(r0+rGSO);
v1Orbit = sqrt(muEarth/r0);
v2Orbit = sqrt(muEarth/rGSO);
v1Transfer = sqrt(2*(muEarth/r0 + eTransfer));
v2Transfer = sqrt(2*(muEarth/rGSO + eTransfer));
tTransfer = sqrt(rTransfer^3/muEarth);

% dV maneuver calcuations
dV1 = v1Transfer - v1Orbit;
dV2 = v2Orbit-v1Transfer;
dVTotal = abs(dV1) + abs(dV2);
fprintf('\nTotal delta V for Hohmann Transfer: %.2f km/s\n',dVTotal/1e3);
fprintf('Time to reach GSO with Hohmann Transfer: %.2f hours\n',tTransfer/3600);

% Plot Hohmann Transfer
figure;
thetaTransfer = linspace(0,pi,500);
xHohmann = rTransfer*cos(thetaTransfer) - (rTransfer-r0);
yHohmann = rTransfer*sin(thetaTransfer);
plot(xInt/rEarth,yInt/rEarth,'g',xFinal/rEarth,yFinal/rEarth,'r');
hold on;
plot(xHohmann/rEarth,yHohmann/rEarth,'color','#0072BD');
grid on;
axis equal;
title('Hohmann Transfer from LEO to GSO');
xlabel('x [Earth Radii]');
ylabel('y [Earth Radii]');

%}
