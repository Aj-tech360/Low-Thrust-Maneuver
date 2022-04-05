function [dR,isterminal,direction] = ltmOdeEventHandler(t,y,r0,rDesired)
    % ltmOdeEventHandler terminates the ode45 function when a desired
    % orbital radius is reached
    %   @param t            -> time vector from ode45 function (not used)
    %   @param y            -> 4x1 array of the states (y = [rho; A; B; theta])
    %   @param r0           -> initial orbit radius in meters
    %   @param rDesired     -> desired orbit radius in meters
    %   @returns dR         -> difference in current and desried orbit radius
    %   @returns isterminal -> stops ode45 integration
    %   @returns direction  -> approaches dR from either side
    % --------------------------------------------------------------------------------
    
    dR = r0*y(1) - rDesired; % stops when this value equals 0
    isterminal = 1;
    direction = 0;
    
end