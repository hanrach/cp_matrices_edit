% explicit solution of heat eq  u_t = u_xx + u_yy, with periodic BCs
% Method-of-lines version
clear;
format long;
errors = [];
gridsize = [1/5, 1/10, 1/15, 1/10/2,];
for hx = gridsize
% Grid and initial data:
% hx = 1/10/4;
 hy = hx;
x1d = -2:hx:(2-hx);
y1d = -2:hy:(2-hy);
%x1d = linspace(0, 1, 21);
%y1d = linspace(0, 1, 41);
%hmin = min(x1d(2)-x1d(1), y1d(2)-y1d(1));
k = .25*(min(hx, hy))^2;    % time step (try .5, what happens?)

[x, y] = meshgrid(x1d, y1d);

% closest point function
R = 1;
[th, r] = cart2pol(x,y);
[cpx, cpy] = pol2cart(th,R);
cpx = cpx(:); cpy = cpy(:);

  theta = linspace(0, 2*pi, 2000);
  xp = R*cos(theta);
  yp = R*sin(theta);
  
% initial condition (must be periodic)

% u0 = sin(2*pi*x) .* sin(2*pi*y);
[th, r] = cart2pol(x,y);
u0 = cos(2*th);

% exact soln, for comparison
% uexact = @(t,x,y) exp(-8*pi^2*t) .* sin(2*pi*x) .* sin(2*pi*y);
uexact = @(t, theta) exp(-4*t)*cos(2*theta);


%Sparse matrix to execute finite difference operation:
[Dxx, Dyy, Dxc, Dyc, Dxb, Dyb, Dxf, Dyf, Dxyc] = ...
  diff2d_matrices(x1d, y1d, 0, 'p');

Tf = 1/8;
% adjust either final time or time-step to have integer steps
numsteps = ceil(Tf / k);
%k = Tf / numsteps
Tf = k*numsteps;

u = u0(:);

% Time-stepping:
% E = interp2(x,y,reshape(u, length(x1d), length(x1d)),cpx,cpy, 'cubic');
for n = 1:numsteps
   unew = u + k*(Dxx*u + Dyy*u);
   unew = interp2(x,y,reshape(unew, length(x1d), length(x1d)),cpx,cpy, 'cubic');
  % cp extensions
   u = unew;
  
  % 4b

  
  unew_theta =  interp2(x,y,reshape(unew, length(x1d), length(x1d)),xp,yp, 'cubic');
  u_theta = unew_theta;
  
  uexact_val = uexact(n*k, theta);
  if (mod(n,50) ==0)
      figure(1);
      plot(theta, u_theta, theta, uexact(n*k,theta));
      drawnow;
  
%   uplot_theta = reshape(u_theta, length(y1d), length(x1d));
%   pcolor(xp,yp,uplot_theta);
  
      figure(2);
      uplot = reshape(u, length(y1d), length(x1d));
      pcolor(x,y,uplot);
      axis equal; axis tight
      shading flat
      %colorbar
      xlabel('x'); ylabel('y');
      drawnow;
      
%       figure(3);
%       uplot_theta = reshape(u_theta, length(xp), length(yp));
%       pcolor(xp,yp,uplot_theta);
%       axis equal; axis tight
%       shading flat
%       %colorbar
%       xlabel('xp'); ylabel('yp');
%       drawnow;
  end
     error = uexact_val - u_theta;
     error_inf = norm(error, inf);
  
  % 4c??
  
end
errors = [errors error_inf];
end

figure(3);
plot(gridsize, errors);
drawnow;

figure(4);
octave_errors = [6.56214468573912e-03 ,  9.33487539523270e-04  , 4.32126862846216e-04, 2.79674066142155e-04];

matlab_errors = vpa(errors);
plot(gridsize, octave_errors, gridsize, matlab_errors);


