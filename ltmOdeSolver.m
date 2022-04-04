function [t,y] = ltmOdeSolver(t,y,r0,v)
    % y = [rho A B theta]

    % Calculate gravitational acceleration 
    g0 = 9.81;
    muEarth = 3.986e14;
    g = g0*(1/y(1))^2;
    
    % Calcuate the derivatives
    dTaudt = sqrt(g0/r0);
    dPdT = y(2)*dTaudt;
    dAdT = dTaudt*(y(3)^2 - y(1))/y(1)^3;
    dBdT = v*y(1)*dTaudt;
    d0dt = sqrt(muEarth/(y(1)*r0)^3);

    % Return the needed dertivatives
    y = [dPdT,dAdT,dBdT,d0dt];
    
end