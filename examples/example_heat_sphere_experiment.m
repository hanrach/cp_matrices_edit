%% Heat equation on a sphere
% cp_matrices is a folder of useful functions to make implementing the
% closest point method easier. These include closest point extension
% matrices, and differentiation matrices.
%
% This example solves the heat equation on the surface of a sphere,
% with initial conditions u = cos(4*theta)

u_err = [];
% diff_plots = [];
prods=[];
base=0.5;
prod =1;
for i = 1:2
    prod=base*prod;
    prods(i) = prod;
end
xdx = 0.1.*prods; xdx =[xdx 0.1 0.15 0.2]; xdx = sort(xdx);
%xdx = (0.1/5):(0.1/5):0.;
for dx = xdx
    dx
    %%
    % 3D example on a sphere
    % Construct a grid in the embedding space

    %dx = 0.1;                   % grid size


    % make vectors of x, y, positions of the grid
    x1d = (-3.0:dx:3.0)';
    y1d = x1d;
    z1d = x1d;

    nx = length(x1d);
    ny = length(y1d);
    nz = length(z1d);

%     [th1d, phi1d, r1d] = cart2sph(x1d, y1d,z1d);
%     [x1dr, y1dr, z1dr] = sph2cart(th1d - angle, phi1d, r1d);


    %% Find closest points on the surface
    % For each point (x,y,z), we store the closest point on the sphere
    % (cpx,cpy,cpz)

    % meshgrid is only needed for finding the closest points, not afterwards
    [xx yy zz] = meshgrid(x1d, y1d, z1d);
    % function cpSphere for finding the closest points on a sphere
    [cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
    % make into vectors
    cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

%     %% Find closest points on the surface - rotated
%     % For each point (x,y,z), we store the closest point on the sphere
%     % (cpx,cpy,cpz)
% 
%     % meshgrid is only needed for finding the closest points, not afterwards
%     [xxr yyr zzr] = meshgrid(x1dr, y1dr, z1dr);
%     % function cpSphere for finding the closest points on a sphere
%     [cpxr, cpyr, cpzr, distr] = cpSphere(xxr,yyr,zzr);
%     % make into vectors
%     cpxgr = cpxr(:); cpygr = cpyr(:); cpzgr = cpzr(:);
    %% Banding: do calculation in a narrow band around the sphere
    dim = 3;    % dimension
    p = 3;      % interpolation degree
    order = 2;  % Laplacian order
    % "band" is a vector of the indices of the points in the computation
    % band.  The formula for bw is found in [Ruuth & Merriman 2008] and
    % the 1.0001 is a safety factor.
    bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
    band = find(abs(dist) <= bw*dx);

    % store closest points in the band;
    cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);
    xg = xx(band); yg = yy(band); zg = zz(band);

%     % store closest points in the band; - rotated
%     cpxgr = cpxgr(band); cpygr = cpygr(band); cpzgr = cpzgr(band);
%     xgr = xxr(band); ygr = yyr(band); zg = zzr(band);
    %% Function u in the embedding space
    % u is a function defined on the grid (eg heat if solving the heat
    % equation)

    % assign some initial value (using initial value of cos (8*theta))
    [th, phi, r] = cart2sph(xx,yy,zz);
    u = cos(3*th + pi/2);

    %assign rotated initial value
    angle = pi/4;
%     [xxr, yyr, zzr] = sph2cart(th, phi + angle, r);
%     [thr,phir, rr] = cart2sph(xxr,yyr,zzr);
    u_rotate = cos(3*(th+angle)+pi/2);

    % this makes u into a vector, containing only points in the band
    u = u(band);

    initialu = u;       % store initial value

    % makes u_rotate into a vector, containing only points in band
    u_rotate = u_rotate(band);
    initialu_rotate = u_rotate;
    %% Construct an interpolation matrix for closest point

    % This creates a matrix which interpolates data from the grid x1d y1d,
    % onto the points cpx cpy.

    disp('Constructing interpolation and laplacian matrices');

    E = interp3_matrix(x1d, y1d, z1d, cpxg, cpyg, cpzg, p, band);

%     Er = interp3_matrix(x1dr, y1dr, z1dr, cpxgr, cpygr, cpzgr, p, band);
    % e.g., closest point extension:
    %u = E*u;

    %% Create Laplacian matrix for heat equation


    L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);
%     Lr = laplacian_3d_matrix(x1dr,y1dr,z1dr, order, band, band);

    %% Construct an interpolation matrix for plotting on sphere

    % plotting grid on sphere, based on parameterization
    paramf = @paramSphere;
    [xp,yp,zp] = sphere(64);
    [tp, pp, rp] = cart2sph(xp, yp, zp);
    [xpr, ypr, zpr] = sph2cart(tp-angle, pp, rp);
    xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);
    [th_plot, phi_plot, r] = cart2sph(xp1,yp1,zp1);
    % Eplot is a matrix which interpolations data onto the plotting grid
    Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

    %%Construct rotated interpolation matrix for plotting on sphere to rotate
    %%back to Eplot
    xp1r = xpr(:); yp1r = ypr(:); zp1r = zpr(:);
    [th_plot_t, phi_plot_t, r_t] = cart2sph(xp1r, yp1r,zp1r);
    Eturn = interp3_matrix(x1d,y1d,z1d, xp1r, yp1r, zp1r,p,band);

    figure(1); clf;

    %% Time-stepping for the heat equation

    Tf = 0.1;
    dt = 0.1*dx^2;
    numtimesteps = ceil(Tf/dt)
    % adjust for integer number of steps
    dt = Tf / numtimesteps
    tic
    for kt = 1:numtimesteps
        % explicit Euler timestepping
        unew = u + dt*(L*u);
        % closest point extension
        u = E*unew;

        unew_rotate = u_rotate +dt*(L*u_rotate);
        u_rotate = E*unew_rotate;
        t = kt*dt;
        if (mod(kt,100) == 0) || (kt < 10) || (kt == numtimesteps)
          %figure(1);
            sphplot = Eplot*u;

            err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
            [t dt dx err]

            sphplot = reshape(sphplot, size(xp));
            surf(xp, yp, zp, sphplot);
            title( ['soln at time ' num2str(t) ', kt= ' num2str(kt)] );
            xlabel('x'); ylabel('y'); zlabel('z');
            %caxis([-1.05 1.05]);   % lock color axis
            axis equal; shading interp;
            %      if ~exist OCTAVE_VERSION camlight left;
            colorbar;
            drawnow(); pause(0);
    %      end
            sphplot_rotate = Eturn*u_rotate;

            err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot_rotate,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
            [t dt dx err]

            sphplot_rotate = reshape(sphplot_rotate, size(xp));
            surf(xp, yp, zp, sphplot_rotate);
            title( ['soln at time ' num2str(t) ', kt= ' num2str(kt) 'rotated'] );
            xlabel('x'); ylabel('y'); zlabel('z');
            %caxis([-1.05 1.05]);   % lock color axis
            axis equal; shading interp;
            %      if ~exist OCTAVE_VERSION camlight left;
            colorbar;
            drawnow(); pause(0);
    %      end
           
        end
    end
    error = norm(sphplot-sphplot_rotate, 2)
    u_err = [u_err error];
    figure(2); clf;
    diff_plot = surf(xp, yp, zp, reshape(sphplot-sphplot_rotate, size(xp)));
    title('difference error plot');
    shading flat;
end

%%plot error of original u and u_rotate turned back to original plotting
%%grid
figure(3); clf;
plot(log(xdx),log(u_err),'.-');
title('error vs. grid size');
xlabel('log h'); ylabel('log error');

