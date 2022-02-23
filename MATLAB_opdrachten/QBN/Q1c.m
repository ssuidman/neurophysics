clear all ;
d = 0.5*10^-3 ;
x = linspace(-2*d,2*d,100) ;
y = linspace(-2*d,2*d,100) ;
[X,Y] = meshgrid(x,y) ;
e0 = 8.85*10^-12 ;

I1 = 1*10^-9 ; 
A1 = I1./(4*pi*e0)./(X.^2+(Y+d).^2).^(3/2) ; 
Ex1 = A1 .* X ; 
Ey1 = A1 .* (Y+d) ; 

I2 = -1*10^-9 ; 
A2 = I2./(4*pi*e0)./(X.^2+(Y-d).^2).^(3/2) ; 
Ex2 = A2 .* X ; 
Ey2 = A2 .* (Y-d) ; 

Ex = Ex1 + Ex2 ;
Ey = Ey1 + Ey2 ;

phi1 = I1./(4*pi*e0.*(X.^2+(Y+d).^2)) ;
phi2 = I2./(4*pi*e0.*(X.^2+(Y-d).^2)) ;
phi = phi1 + phi2 ;

figure(1) ;
quiver(X,Y,tanh(Ex/10^8),tanh(Ey/10^8)) ; 
hold on
contour(X,Y,tanh(phi/10^8)) ;
