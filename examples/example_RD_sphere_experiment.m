%% Reaction-diffusion on a sphere

cpf = @cpSphere;
paramf = @paramSphere;
%cpf = @cpEllipsoid;
%paramf = @paramEllipsoid;

loaddata = 1;

errors_rotate_back = [];
errors_rotate_org = [];
errors_comp_grid = [];
errors_comp_grid_rotate =[];

%%initialize grid sizes to test
h=0.005
xdx = [0.025, 0.05, 0.1, 0.15]
xdx = 0.025:h:0.1;
for dx = xdx
    if (loaddata == 1)
      % grid size
      dx_original = 0.05;
      angle = pi/4
           
      % make vectors of x, y, z positions of the grid
      x1d = (-2.0:dx:2.0)';
      y1d = x1d;
      z1d = x1d;
      nx = length(x1d);
      ny = length(y1d);
      nz = length(z1d);

      % meshgrid is only needed for finding the closest points, not afterwards
      [xx yy zz] = meshgrid(x1d, y1d, z1d);

      [cpx, cpy, cpz, dist] = cpf(xx,yy,zz);
      cpx = cpx(:); cpy = cpy(:); cpz = cpz(:);

      %% Banding: do calculation in a narrow band around the surface
      dim = 3;  % dimension
      p = 3;    % interpolation order
      % "band" is a vector of the indices of the points in the computation
      % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
      % the 1.0001 is a safety factor.
      bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((1+(p+1)/2)^2));
      band = find(abs(dist) <= bw*dx);

      % store closest points in the band;
      cpx = cpx(band); cpy = cpy(band); cpz = cpz(band);
      x = xx(band); y = yy(band); z = zz(band);


      %% discrete operators
      disp('building laplacian and interp matrices');
      L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band,band);
      E = interp3_matrix(x1d,y1d,z1d, cpx, cpy, cpz, p, band);
      I = speye(size(E));

      %% plotting grid
      [xp,yp,zp] = paramf(256);

      % Eplot is a matrix which interpolations data onto the plotting grid
      Eplot = interp3_matrix(x1d, y1d, z1d, xp(:), yp(:), zp(:), p, band);

      %% build spherical coordinate system
    %   r = sqrt(x.^2+y.^2+z.^2);
    %   theta = acos(z./r);
    %   phi =atan(y./x);

    [theta, phi, r]= cart2sph(x,y,z);
    [x2, y2, z2] = sph2cart(theta - angle, phi, r);
    
    [thetap,phip,rp] = cart2sph(xp,yp,zp);
    [x2p, y2p, z2p] = sph2cart(thetap - angle, phip, rp);

    Eturn_plot =interp3_matrix(x1d, y1d, z1d, x2p(:), y2p(:), z2p(:), p, band);
    E_comp = interp3_matrix(x1d, y1d, z1d, x, y, z, p, band);
    Eturn = interp3_matrix(x1d, y1d, z1d, x2, y2, z2, p, band);
    end

    figure(1); clf;

    % u_t = f(u,g) + nuu*Lap u
    % v_t = g(u,g) + nuv*Lap u

    % parameters and functions for Gray--Scott
    FF = 0.054;  kk = 0.063;  nuu = 1/3600;  nuv = nuu/3;
    f = @(u,v) (-u.*v.*v  +  FF*(1-u));
    g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);
    


    %% initial conditions1 - small perturbation from steady state

    xrand = rand(size(x));
    yrand = rand(size(y));
    zrand = rand(size(z));

    [theta_r, phi_r, r_r] = cart2sph(xrand,yrand,zrand);
    [xrandt, yrandt, zrandt] = sph2cart(theta_r+ pi/4, phi_r, r_r);

    pert = (1/4)*exp(-(10*(z-.1)).^2) + (1/4)*sin(3*(theta)); %(1/4)*sin(7*(theta)); %+ 0.5*xrand;
    u0 = 1-pert;  v0 = 0.5*pert;
    u = u0;  v = v0;
    
    pert_rotate = (1/4)*exp(-(10*(z-.1)).^2) + (1/4)*sin(3*(theta+angle)); 
    u0_rotate = 1-pert_rotate; v0_rotate = 0.5*pert_rotate;
    u_rotate = u0_rotate; v_rotate = v0_rotate;

    Tf = 100;
    dt = .2 * (1/max(nuu,nuv)) * dx^2;
    numtimesteps = ceil(Tf/dt);
    % adjust for integer number of steps
    dt = Tf / numtimesteps;


    figure(1);
    sphplot = Eplot*u;
    sphplot = reshape(sphplot, size(xp));
    Hplot = surf(xp, yp, zp, sphplot);
    title('initial u')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    view(-10, 60)
    axis off;
    shading interp
    camlight left
    colorbar

    sphplot_rotate = Eturn_plot*u_rotate;
    sphplot_rotate = reshape(sphplot_rotate, size(x2p));
    Hplot_rotate = surf(x2p, y2p, z2p, sphplot_rotate);
    title('initial u rotate')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    view(-10, 60)
    axis off;
    shading interp
    camlight left
    colorbar
    %% Method-of-lines approach
    % See [vonGlehn/Macdonald/Maerz 2013]
    %lambda = 6*max(nuu,nuv)/(dx^2);
    %Au = nuu*(E*L) - lambda*(I-E);
    %Av = nuv*(E*L) - lambda*(I-E);


    for kt = 1:numtimesteps
      %% MOL: explicit Euler timestepping
      %unew = u + dt*( E*f(u,v) + Au*u );
      %vnew = v + dt*( E*g(u,v) + Av*v );
      %u = unew;
      %v = vnew;
      %% MOL: without precomputing matrices
      %rhsu = nuu*(L*u) + f(u,v);
      %rhsv = nuv*(L*v) + g(u,v);
      %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
      %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
      %u = unew;
      %v = vnew;

      %% Ruuth-Merriman
      rhsu = nuu*(L*u) + f(u,v);
      rhsv = nuv*(L*v) + g(u,v);
      unew = u + dt*rhsu;
      vnew = v + dt*rhsv;
      u = E*unew;
      v = E*vnew;
      
      rhsu_r = nuu*(L*u_rotate) + f(u_rotate,v_rotate);
      rhsv_r = nuv*(L*v_rotate) + g(u_rotate,v_rotate);
      unew_r = u_rotate + dt*rhsu_r;
      vnew_r = v_rotate + dt*rhsv_r;
      u_rotate = E*unew_r;
      v_rotate = E*vnew_r;

      t = kt*dt;

      if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
        disp([kt t]);
        sphplot = Eplot*u;

        sphplot = reshape(sphplot, size(xp));
        surf(xp, yp, zp, sphplot);
%         set(0, 'CurrentFigure', 1);
%         set(Hplot, 'CData', sphplot);
        title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
        xlabel('x'); ylabel('y'); zlabel('z');
        axis equal
        view(-10, 60)
        axis off;
        shading interp
        camlight left
        colorbar
        drawnow(); 
        
        
        disp([kt t]);
        sphplot_rotate = Eturn_plot*u_rotate;

        sphplot_rotate = reshape(sphplot_rotate, size(x2p));
%         set(0, 'CurrentFigure', 2);
%         set(Hplot_rotate, 'C2Data', sphplot_rotate);
        surf(xp, yp, zp, sphplot_rotate)
        title( ['u at time ' num2str(t) ', kt= ' num2str(kt) ', rotated'] );
        xlabel('x'); ylabel('y'); zlabel('z');
        axis equal
        view(-10, 60)
        axis off;
        shading interp
        camlight left
        colorbar
        drawnow();

      end
    end
        
    %% u_orgiginal - u_rotate rotated -pi/4 on computation grid
%     u_error = norm(u_original-Eturn*u_rotate,inf);
    %%u_original*Eturn - u_rotate*Eplot\
    
    %% 
    u_error_rotate_back = norm(Eplot*u - Eturn_plot*u_rotate, inf);
    u_error_rotate_org = norm(Eturn_plot*u - Eplot*u_rotate, inf);
    u_error_comp_grid = norm(E_comp*u - Eturn*u_rotate, inf);
    u_error_comp_grid_rotate = norm(u - Eturn*u_rotate, inf);
    errors_rotate_back = [errors_rotate_back u_error_rotate_back];
    errors_rotate_org = [errors_rotate_org u_error_rotate_org];
    errors_comp_grid = [errors_comp_grid u_error_comp_grid];
    errors_comp_grid_rotate = [errors_comp_grid_rotate u_error_comp_grid_rotate];
end

% 
%     figure(2); clf;
% 
%     % u_t = f(u,g) + nuu*Lap u
%     % v_t = g(u,g) + nuv*Lap u
% 
%     % parameters and functions for Gray--Scott
%     FF = 0.054;  kk = 0.063;  nuu = 1/3600;  nuv = nuu/3;
%     f = @(u,v) (-u.*v.*v  +  FF*(1-u));
%     g = @(u,v) ( u.*v.*v  -  (FF+kk)*v);
%     %% initial conditions2 - small perturbation from steady state
% 
%     pert = (1/4)*exp(-(10*(z-.1)).^2) + (1/4)*exp(-(10*(x2-.1)).^2); %(1/4)*sin(7*(theta+pi/4));%+0.5*xrandt;
%     u0 = 1-pert;  v0 = 0.5*pert;
%     u = u0;  v = v0;
% 
%     Tf =30;
%     dt = .2 * (1/max(nuu,nuv)) * dx^2;
%     numtimesteps = ceil(Tf/dt);
%     % adjust for integer number of steps
%     dt = Tf / numtimesteps;
% 
%     figure(2);
%     sphplot = Eplot*u;
% 
%     sphplot = reshape(sphplot, size(x2p));
%     Hplot = surf(x2p, y2p, z2p, sphplot);
%     title('initial u')
%     xlabel('x'); ylabel('y'); zlabel('z');
%     axis equal
%     view(-10, 60)
%     axis off;
%     shading interp
%     camlight left
%     colorbar
% 
%     %% Method-of-lines approach
%     % See [vonGlehn/Macdonald/Maerz 2013]
%     %lambda = 6*max(nuu,nuv)/(dx^2);
%     %Au = nuu*(E*L) - lambda*(I-E);
%     %Av = nuv*(E*L) - lambda*(I-E);
% 
%     for kt = 1:numtimesteps
%       %% MOL: explicit Euler timestepping
%       %unew = u + dt*( E*f(u,v) + Au*u );
%       %vnew = v + dt*( E*g(u,v) + Av*v );
%       %u = unew;
%       %v = vnew;
%       %% MOL: without precomputing matrices
%       %rhsu = nuu*(L*u) + f(u,v);
%       %rhsv = nuv*(L*v) + g(u,v);
%       %unew = u + dt*( E*rhsu - lambda*(u - E*u) );
%       %vnew = v + dt*( E*rhsv - lambda*(v - E*v) );
%       %u = unew;
%       %v = vnew;
% 
%       %% Ruuth-Merriman
%       rhsu = nuu*(L*u) + f(u,v);
%       rhsv = nuv*(L*v) + g(u,v);
%       unew = u + dt*rhsu;
%       vnew = v + dt*rhsv;
%       u = E*unew;
%       v = E*vnew;
% 
%       t = kt*dt;
% 
%       if ( (mod(kt,25)==0) || (kt<=10) || (kt==numtimesteps) )
%         disp([kt t]);
%         sphplot = Eplot*u;
%         sphplot = reshape(sphplot, size(x2p));
%         set(0, 'CurrentFigure', 2);
%         set(Hplot, 'CData', sphplot);
%         title( ['u at time ' num2str(t) ', kt= ' num2str(kt)] );
%         drawnow;
%       end
%     end
%     u_rotate = u;
%     
%     %% u_orgiginal - u_rotate rotated -pi/4 on computation grid
%     u_error = norm(u_original-Eturn*u_rotate,inf);
%     %%u_original*Eturn - u_rotate*Eplot\
%     
%     %% u_original rotated pi/4 - u_original both on plotting grid
%     u_error = norm(Eturn_plot*u_original - Eplot*u_rotate, inf);
%     u_err = [u_err, u_error]
%     figure(4);
%     
% end


%%plot grid size vs. u_err
figure(3);
delta_x_plot = plot(xdx, errors_rotate_back, '.-');
title('u_rotate on rotated grid vs. grid size(h)');
xlabel('h'); ylabel('norm(Eplot*u - Eturn_plot*u_rotate, inf)');

figure(4);
delta_x_plot = plot(xdx, errors_comp_grid, '.-');
title('u_rotate on rotated computation and u on comp grid vs. grid size(h)');
xlabel('h'); ylabel('norm(E_comp*u - Eturn*u_rotate, inf)');

figure(5);
delta_x_plot = plot(xdx, errors_comp_grid_rotate, '.-');
title('u_rotate on rotated computation vs. grid size(h)');
xlabel('h'); ylabel('norm(u - Eturn*u_rotate, inf)');

%%plot the difference - only captures the last iteration of the for loop
figure(6);
difference_plot = surf(xp,yp,zp,reshape(Eturn_plot*u-Eplot*u_rotate,size(xp))); 
shading flat


