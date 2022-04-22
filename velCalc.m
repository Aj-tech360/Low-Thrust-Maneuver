function vel = velCalc(y,r0,v0,g0)
    % velCal solves for the dimensional velocity for a given LTM orbit 
    % velCal returns the dimensional velocity in m/s
    %   @param y       -> 4x1 array of the states (y = [rho A B theta])
    %   @param r0      -> initial orbit radius in meters
    %   @param v0      -> initial orbital velcity in m/s
    %   @returns vel   -> dimensional velocity in m/s 
    % --------------------------------------------------------------------------------

    % Values from state array
    rho = y(:,1);
    A = y(:,2);
    B = y(:,3);

    % Dimensional velocity calculation
    r = r0*rho;
    rDot = A*sqrt(r0*g0);
    thetaDot = (B./(rho.^2))*sqrt(g0/r0);
    vRadSq = rDot.^2;
    vTanSq = (r.*thetaDot).^2;
    velDim = sqrt(vRadSq+vTanSq);

    % Return spacecraft velocity
    vel =  velDim;
end







