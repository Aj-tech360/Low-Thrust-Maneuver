function dydt = ltmOdeSolver(~,y,r0,g0,v)
    % ltmOdeSolver solves for the first-order ODEs for use in the low thrust maneuver 
    % ODE for electrical rocekt engines
    % ltmOdeSovler returns the first-order ODEs for use in the ode45 function 
    %   @param y       -> 4x1 array of the states (y = [rho A B theta])
    %   @param r0      -> initial orbit radius in meters
    %   @param g0      -> gravitational acceleration at initial orbit radius in
    %                     meters/second^2
    %   @param v       -> dimensionless thrust factor
    %   @returns dydt  -> 4x1 array of the state derivatives (dydt) 
    % --------------------------------------------------------------------------------

    % Calcuate the derivatives
    dTaudt = sqrt(g0/r0);
    dPdt = y(2)*dTaudt;                            
    dAdt = dTaudt*((y(3)^2 - y(1))/(y(1)^3));    
    dBdt = v*y(1)*dTaudt;       
    dAdT = dAdt*sqrt(r0/g0);
    d0dT = sqrt((1/y(1))*(dAdT+(1/(y(1))^2)));
    d0dt = d0dT*dTaudt;
    
    % Return the needed dertivatives
    dydt = [dPdt;dAdt;dBdt;d0dt];
end