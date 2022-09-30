%%
%1.4)
x = linspace(0,0.1,1000)
w_car = 1000*2*pi
w_mod = 100*2*pi
y1 = (1+sin(w_mod*x)).*sin(w_car*x)
y2 = max(0,y1)

freq = linspace(0,5000,length(x)/2+1)
ft1 = abs(fft(y1,1024))
ft1 = ft1(1:length(x)/2+1)
ft2 = abs(fft(y2,1024))
ft2 = ft2(1:length(x)/2+1)
figure(1)
plot(freq,ft1,'color','b')
hold on
plot(freq,ft2,'color','r')
legend('normal','rectified')
xlim([0,5000])
xlabel('Frequency (Hz)')
title('Frequency spectra of normal and rectified modulated signal')

%It can be seen that there are multiple spikes in the graph. The spike at
%f_car=1000Hz, but also f_mod=100Hz, f_car+f_mod=1100Hz,
%f_car-f_mod=900Hz can be seen. This is the same as in 6.17. For the
%rectified signal there are also peaks at higher frequencies, such as
%2*f_car+/-f_mod and 4*f_car+/-f_mod


