%%rotation invariance rotate a vector
dx = 0.5;
x1d = (-10:dx:10);
y1d = x1d;
z1d = x1d;

[x y] = meshgrid(x1d, y1d);

%%simpler testing



%% build spherical coordinate system
[theta,r] = cart2pol(x,y);

z0 = x.^2+y.^2;
z1 = 4*(r.*sin(theta+pi/4));%-25+exp(10*sin(r.*cos(theta+pi/4))); 
z2 = 4*(r.*sin(theta+pi/2));%-25+exp(10*sin(r.*cos(theta+pi/2)));
z3 = 4*(r.*sin(theta));%-25+exp(10*sin(r.*cos(theta+pi)));
%z=4*sin(y)-25+exp(10*sin(x));


figure;
subplot(4,1,1) 
%surf(x,y,z);
surf(r,theta,z1);


subplot(4,1,2) 
%surf(x,y,z);
surf(r,theta,z2);

subplot(4,1,3) 
%surf(x,y,z);
surf(z0,x,y);

subplot(4,1,4)
h0 = surf(z0,x,y);
direction = [0 0 1];
rotate(h0,direction,pi);

% R = [cos(2*pi), -sin(2*pi); sin(2*pi), cos(2*pi)];
% 
% frot = R*f;


