function dydt = ltmOdeSolver(t,y,r0,v)
    % ltmOdeSolver solves for the first-order ODEs for use in the low thrust maneuver 
    % ODE for electrical rocekt engines
    % ltmOdeSovler returns the first-order ODEs for use in the ode45 function 
    %   @param t   -> time vector from ode45 function (not used)
    %   @param y   -> 4x1 array of the states (y = [rho A B theta])
    %   @param r0  -> initial orbit radius in meters
    %   @param v   -> 
    %   @returns   -> 4x1 array of the state derivatives (dydt) 
    % --------------------------------------------------------------------------------

    % Constants 
    g0 = 9.81;
    muEarth = 3.986e14;

    % Calcuate the derivatives
    dTaudt = sqrt(g0/r0);
    dPdt = y(2)*dTaudt;
    dAdt = dTaudt*((y(3)^2 - y(1))/y(1)^3);
    dBdt = v*y(1)*dTaudt;
    d0dt = sqrt(muEarth/(y(1)*r0)^3);

    % Return the needed dertivatives
    dydt = [dPdt;dAdt;dBdt;d0dt];
end