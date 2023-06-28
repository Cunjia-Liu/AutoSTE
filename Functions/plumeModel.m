function C = plumeModel(s, p)
% Plume dispersino model
% s is the source term (either a structure or array)
% p contains the locations of the concentration readings
% C is the concentration level from the plume model

% Dispersion model assumes isotropic diffusion from the source, taken from 
% Vergassola et al. (2007) and also in Hutchinson et al. (2018)
% see IP model in https://onlinelibrary.wiley.com/doi/epdf/10.1002/rob.21844


if isstruct(s)
    s.D=s.ci;
    s.t = s.cii;
    lamda = sqrt((s.D.*s.t)./(1+ (s.u.^2.*s.t)./(4*s.D)));
    
    module_dist = sqrt(((s.x-p.x_matrix)).^2 + ((s.y-p.y_matrix)).^2 + ((s.z-p.z_matrix)).^2);

    module_dist(module_dist<1e-5)=1e-5;
    
    C = s.Q./(4*pi.*s.D.*module_dist).*exp((-(p.x_matrix-s.x).*s.u.*cos(s.phi)./(2.*s.D)) + (-(p.y_matrix-s.y).*s.u.*sin(s.phi)./(2.*s.D)) + (-1*module_dist./lamda));
else
    x = s(1,:)';
    y = s(2,:)';
    z = s(3,:)';
    Q = s(4,:)';
    u = s(5,:)';
    phi = s(6,:)';
    D = s(7,:)';
    t = s(8,:)';

    lamda = sqrt((D.*t)./(1+ (u.^2.*t)./(4*D)));
    
    module_dist = sqrt(((x-p.x_matrix)).^2 + ((y-p.y_matrix)).^2 + ((z-p.z_matrix)).^2);

    module_dist(module_dist<1e-5)=1e-5;
    
    C = Q./(4*pi.*D.*module_dist).*exp((-(p.x_matrix-x).*u.*cos(phi)./(2.*D)) + (-(p.y_matrix-y).*u.*sin(phi)./(2.*D)) + (-1*module_dist./lamda));

end





