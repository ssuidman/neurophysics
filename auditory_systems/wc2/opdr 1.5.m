%%
%1.5)
%f_c=1000Hz, f_m=100Hz, b=1

t1 = linspace(0,0.1,1000)
w1 = linspace(0,5000,length(t1)/2+1)

w_c1 = 1000*2*pi
w_m1 = 100*2*pi
b1 = 1
s1 = sin(w_c1*t1-b1*cos(w_m1*t1))

fts1 = abs(fft(s1,1024))
fts1 = fts1(1:length(t1)/2+1)

figure(2)
plot(w1,fts1)
xlabel('Frequency (Hz)')
title('Frequency spectra, f_{carrier}=1000Hz, f_{modulated}=100Hz,b=1')


%%
%f_c=3000Hz, f_m=200Hz, b=3

t2 = linspace(0,0.1,1000)
w2 = linspace(0,5000,length(t2)/2+1)

w_c2 = 3000*2*pi
w_m2 = 200*2*pi
b2 = 3
s2 = sin(w_c2*t2-b2*cos(w_m2*t2))

fts2 = abs(fft(s2,1024))
fts2 = fts2(1:length(t2)/2+1)

figure(3)
plot(w2,fts2)
xlabel('Frequency (Hz)')
title('Frequency spectra, f_{carrier}=3000Hz, f_{modulated}=200Hz,b=3')


%This signal gives for b=0 also only a peak at f=f_car=1000Hz (which is 
%also the middle peak if b is not 0), however if you
%increase b, more and more peaks arise to the sides of it with +/-f_mod,
%+/-2*f_mod, +/-3*f_mod, etc. So if f_mod becomes bigger then the spikes of
%the side peaks become more apart of eachother. 












