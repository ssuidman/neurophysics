clear all ; 
d = 0.5 ; 
% d = 0.5*10^-3 ; 
x = linspace(-2*d,2*d,100) ; 
y = linspace(-2*d,2*d,100) ; 
% z = linspace(-2*d,2*d,100) ; 
[X,Y] = meshgrid(x,y) ; 
e0 = 8.85*10^-12 ; 
I = 1*10^-3 ; 

r = sqrt(X.^2+Y.^2) ;
pz = I*(2*d) ;
% Ez = (3*pz*Z.^2 - pz*r.^2) ./ r.^5 ; 
% Ey = (3*pz*Z.*Y) ./ r.^5 ; 
Ey = (3*pz*Y.^2 - pz*r.^2) ./ r.^5 ;
Ex = (3*pz*Y.*X) ./ r.^5 ;

phi = (pz*Y) ./ (4*pi*e0*r.^3) ;

figure(1) ;
contour(X,Y,phi) ;
hold on ;
% quiver(X,Y,Ex,Ey) ; 
quiver(X,Y,Ex./abs(Ex),Ey./abs(Ey)) ; 
