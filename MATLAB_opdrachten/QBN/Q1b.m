clear all ; 
x = linspace(-1,1,100) ; 
y = linspace(-1,1,100) ; 
[X,Y] = meshgrid(x,y) ; 
e0 = 8.85*10^-12 ;
 
I = 1*10^-9 ; 
A = I./(4*pi*e0)./(X.^2+Y.^2).^(3/2) ; 
Ex = A .* X ; 
Ey = A .* Y ; 


phi = I./(4*pi*e0.*(X.^2+Y.^2))

figure(1) ; 
quiver(X,Y,Ex./abs(Ex),Ey./abs(Ey)) ; 

hold on ;

contour(X,Y,phi) ;


