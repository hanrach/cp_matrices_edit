    r_rand = rand(size(r));
    theta_rand = rand(size(theta));
    phi_rand = rand(size(phi));
    r_pert = r_rand.*sin(theta_rand).*cos(phi_rand);
    % pert = 0.5*exp(-(10*(z-.1)).^2); %+ 0.5*randx %+ 0.5*rand(size(y));
    pert_z = 0.5*exp(-(10*((r.*cos(theta))-.1)).^2);%+0.5*(r_pert);
    pert_x = -(10*(r.*sin(theta).*cos(phi))-.1);
    pert =  0.5*exp(-(10*z-.1)).^2;
    x_coord = (r.*sin(theta)).*(cos(phi));   

r_pert = r_rand.*sin(theta_rand+(pi)).*cos(phi_rand);
    % pert = 0.5*exp(-(10*(z-.1)).^2); %+ 0.5*randx %+ 0.5*rand(size(y));
    pert_z = 0.5*exp(-(10*(r.*cos(theta)-.1)).^2);%+0.5*(r_pert);
    pert_x = 0.5*exp(-(10*(r.*sin(theta+(pi/2)).*cos(phi+(pi/2)))-.1)).^2;
    pert = pert_z;
    pert = 0.5*exp(-(10*z2-.1).^2);
    

% R = [cos(2*pi), -sin(2*pi), 0; sin(2*pi), cos(2*pi), 0;0 ,0 ,1];
    % 
    % u0 = 1-pert;  v0 = 0.5*pert; %z0 = z./z;
    % 
    % 
    % for i=1:length(u0)
    % 
    % v_rot(:,i) = [u0(i) v0(i) 1]*R;
    % 
    % end
    % 
    % % pick out the vectors of rotated x- and y-data
    % 
    % 
    % 
    % 
    % u_rotated = v_rot(1,:).';
    % v_rotated = v_rot(2,:).';