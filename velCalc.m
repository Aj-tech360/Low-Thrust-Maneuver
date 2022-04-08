function vel = velCalc(y,r0,v0)
    g0 = 9.81;
    % From state vector
    dtdTau = sqrt(g0/r0)
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