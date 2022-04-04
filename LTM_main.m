% AE 414 Low Thrust Maneuver Project
close all;clear;clc;

% Constants
g0 = 9.81;
rEarth = 6378e3;
muEarth = 3.986e14;

% Given spacecraft/orbit data
rOrbit = 6698e3;
v = 500e-5;
t0 = 0;
tBurn = 172800; % 2 days => s
vOrbit0 = sqrt(muEarth/rOrbit);

% Graviational accleration function 
g = @(r) g0*(r0/r)^2;

% ODE initial conditions
IC = [0 1 1]; %dp/dt p do/dt
nPts = 5000;
tSpan = linspace(t0,tBurn,nPts);

% Call ode45 function 
tau = tSpan*sqrt(g0/rOrbit);

[t,y] = ode45(ltmOdeSolver,tSpan,IC);
