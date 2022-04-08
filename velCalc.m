function vel = velCalc(y,r0,v0)
    % velCal solves for the normalized and dimensional velocity for a given LTM orbit 
    % velCal returns the normalized and dimensional velocity in m/s
    %   @param y       -> 4x1 array of the states (y = [rho A B theta])
    %   @param r0      -> initial orbit radius in meters
    %   @param v0      -> initial orbital velcity in m/s
    %   @returns vel   -> 1x2 array of the normalized velocity and dimensional velocity 
    % --------------------------------------------------------------------------------

    g0 = 9.81;

    % From state vector
    dtdTau = sqrt(g0/r0);
    rho = y(:,1);
    A = y(:,2)*dtdTau;
    %B = y(:,3)*dtdTau;
    %theta = y(:,4);
    %r = r0*rho;

    % Calculate Derivatives
    %rDot = r0*A;
    thetaDot = sqrt((1./rho).*(A + (1./rho.^2)));

    % Return dimensional velocity      %TODO: also return normalized velocity 
    vel = v0+sqrt(r0*g0)*sqrt(A.^2 + rho.^2.*thetaDot.^2);
    


    %sqrt(rDot.^2 + thetaDotSquared.*(r.^2));
    %rDot = r0*y(:,2);
    %thetaDotSquared = ((1./y(:,1)).*(y(:,2)*sqrt(g0/r0) + (1./y(:,1).^2)))*(r0/g0);
end