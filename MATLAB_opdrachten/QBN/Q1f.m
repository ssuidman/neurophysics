clear all ; 
d = 1*10^-3 ; 
x = linspace(-d,d,100) ; 
y = linspace(-d,d,100) ;  
[X,Y] = meshgrid(x,y) ; 
e0 = 8.85*10^-12 ; 
I = 3*10^-9 ; 

r = sqrt(X.^2+Y.^2) ; 
pz = I*(2*d) ; 
Ey = (3*pz*Y.^2 - pz*r.^2) ./ r.^5 ; 
Ex = (3*pz*Y.*X) ./ r.^5 ; 

phi = (pz*Y) ./ (4*pi*e0*r.^3) ; 
r = 1*10^-3
phi_r = phi(100,1) % to calculate the potential at r we must look in our frame to the most right value

disp('potential at r=1mm is:') ;
disp(phi_r) ;

figure(1) ; 
quiver(X,Y,tanh(Ex*10),tanh(Ey*10)) ; 
hold on ;
contour(X,Y,tanh(phi/10^5)) ;
hold on ;
scatter(r,0,100,'filled') ;



