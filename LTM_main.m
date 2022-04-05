% AE 414 Low Thrust Maneuver Project
close all;clear;clc;
 
% Constants
g0 = 9.81;
rEarth = 6378e3;
muEarth = 3.986e14;

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
x = r0*y(:,1).*cos(y(:,4));
y = r0*y(:,1).*sin(y(:,4));

figure;
plot(x/10^3,y/10^3);
grid on;
axis equal;
title('Spacecraft Orbit Over Two Days');
xlabel('x [km]');
ylabel('y [km]');



