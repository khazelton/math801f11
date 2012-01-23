function [x, y, z] = torus (r, n, a)
% TORUS Generate a torus.
% torus (r, n, a) generates a plot of a
% torus with central radius a and
% lateral radius r.
% n controls the number of facets
% on the surface.
% These input variables are optional
% with defaults r = 0.5, n = 30, a = 1.
%
% [x, y, z] = torus(r, n, a) generates
% three (n + 1)-by-(n + 1) matrices so
% that surf (x, y, z) will produce the
% torus.
%
% See also SPHERE, CYLINDER
%
% MATLAB Primer, 6th Edition
% Kermit Sigmon and Timothy A. Davis
% Section 11.5, page 65.

if nargin < 3, a = 1 ; end
if nargin < 2, n = 200 ; end
if nargin < 1, r = 0.4 ; end
theta = pi * (0:2:2*n)/n ;
phi = 2*pi* (0:2:n)'/n ;
[xx1,yy1] = meshgrid(theta,phi);
xx = (a + r*cos(phi)) * cos(theta) ;
yy = (a + r*cos(phi)) * sin(theta) ;
zz = r * sin(phi) * ones(size(theta)) ;
size(xx)
if nargout == 0
    figure;
    mesh (xx, yy, zz,'marker','.','edgecolor','none','facecolor','none','markeredgecolor','b','markersize',2) ;
    ar = (a + r)/sqrt(2) ;
    axis([-ar, ar, -ar, ar, -ar, ar]) ;
    hold on;
    surf(xx(1:20,130:150),yy(1:20,130:150),zz(1:20,30:50),'edgecolor','none');
    surf(xx(91:end,130:150),yy(91:end,130:150),zz(91:end,30:50),'edgecolor','none');
    colormap gray;
    axis equal;
    view(18,34)
else
    x = xx ;
    y = yy ;
    z = zz ;
end
figure;
plot(xx1(:),yy1(:),'b.');
xx1 = xx1(91:end,130:150);
yy1 = yy1(91:end,130:150);
hold on;
plot(xx1(:),yy1(:),'k.');

