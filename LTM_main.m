% Low Thrust Maneuver Project
% Ronak Amin, Benjamin Sites, Christopher Rappole 
% AE 414 – 01 
% Prof. Laksh Narayanaswami 
% April 22, 2022 

close all;clear;clc;
 
% Constants
gEarth = 9.81;
rEarth = 6378e3;
muEarth = 3.986e14;

%% LTM for 2 Days
fprintf('\n\t\tProblem 1\n-------------------------------------------\n');

% Gravitational accelration function
g = @(r) gEarth*(rEarth/r)^2;

% Given spacecraft/orbit data
r0 = 6698e3;
g0 = g(r0);
hGSO = 35786e3;
rGSO = hGSO + rEarth;
v = 500e-5;
t0 = 0;
tBurn = 172800; % 2 days => s
vOrbit0 = sqrt(muEarth/r0);

% ODE initial conditions
IC = [1;0;1;0]; % [rho0 A0 B0 theta0]
nPts = 10000;
tSpan = linspace(t0,tBurn,nPts);

% ode45 function to find rho
[t,y] = ode45(@(t,y) ltmOdeSolver(t,y,r0,g0,v),tSpan,IC);   % y = [rho; A; B; theta]

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
uDim = velCalc(y,r0,vOrbit0,g0);

% Find minimum velocity and dimensional time (in hours)
[minU,idxMinU] = min(uDim);
minUTau = tau(idxMinU);
minUTime = (minUTau*sqrt(r0/g0))/3600; % hours
fprintf('The minimum velocity is %.2f km/s at %.2f hr\n',minU/1e3,minUTime);

% Find velocity at GSO
[~,idxRad] = min(abs(r-rGSO));
uAtGSO = uDim(idxRad);
fprintf('The velocity at GSO altitude is %.2f km/s\n',uAtGSO/1e3);


% Plot velocity vs normalized time
figure;
plot(tau,uDim/1e3);
grid on;
title('Dimensional Velocity of Spacecraft Over 2 Days');
xlabel('Normalized Time');
ylabel('Velocity [km/s]');


%% Spacecraft LTM Orbit Transfer
fprintf('\n\t\tProblem 2\n-------------------------------------------\n');

% Given spacecraft/orbit parameters
v = 2.7e-5;
tSpan = linspace(0,3e7,nPts*10);

% Calculate rho until r0*rho = rGSO
opts = odeset('Events',@(t,y) ltmOdeEventHandler(t,y,r0,rGSO));
[t,y,te,ye,ie] = ode45(@(t,y) ltmOdeSolver(t,y,r0,g0,v),tSpan,IC,opts); % y = [rho; A; B; theta]
transferTime = te/86400; % time to reach orbit in days
fprintf('Time to reach GSO: %.2f days\n',transferTime);

% Calculate delta V 
transferVel = velCalc(y,r0,vOrbit0,g0);
fprintf('The spacecraft''s velocity at GSO altitude is: %.2f km/s\n',transferVel(end)/1e3);
accTransfer = v*g(r0);
dvLtmTransfer = accTransfer*te;
fprintf('Total delta V for LTM Transfer: %.2f km/s\n',dvLtmTransfer/1e3);

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
aTransfer = (r0+rGSO)/2;
eTransfer = -muEarth/(r0+rGSO);
v1Orbit = sqrt(muEarth/r0);
v2Orbit = sqrt(muEarth/rGSO);
v1Transfer = sqrt(2*((muEarth/r0) + eTransfer));
v2Transfer = sqrt(2*((muEarth/rGSO) + eTransfer));
eccTransfer = 1-(r0/aTransfer);
tTransfer = pi*(1/sqrt(muEarth))*aTransfer^(3/2);

% dV maneuver calcuations
dV1 = v1Transfer - v1Orbit;
dV2 = v2Orbit-v2Transfer;
dVTotal = abs(dV1) + abs(dV2);
fprintf('\nTotal delta V for Hohmann Transfer: %.2f km/s\n',dVTotal/1e3);
fprintf('Time to reach GSO with Hohmann Transfer: %.2f hours\n',tTransfer/3600);

% Plot Hohmann Transfer
figure;
thetaTransfer = linspace(0,pi,500);
xHohmann = aTransfer*(cos(thetaTransfer)-eccTransfer);
yHohmann = aTransfer*sqrt(1-eccTransfer^2)*sin(thetaTransfer);
plot(xInt/rEarth,yInt/rEarth,'g',xFinal/rEarth,yFinal/rEarth,'r');
hold on;
plot(xHohmann/rEarth,yHohmann/rEarth,'color','#0072BD');
grid on;
axis equal;
title('Hohmann Transfer from LEO to GSO');
xlabel('x [Earth Radii]');
ylabel('y [Earth Radii]');


