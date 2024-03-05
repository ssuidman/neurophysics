%%
% 1.3a
ft1 = fft(y1,1024)
ft1 = ft1(1:512)
fa1 = abs(ft1)
freq = linspace(0,5000,512)

ft2 = fft(y2,1024)
ft2 = ft2(1:512)
fa2 = abs(ft2)
freq = linspace(0,5000,512)

ft3 = fft(y3,1024)
ft3 = ft3(1:512)
fa3 = abs(ft3)
freq = linspace(0,5000,512)

figure(3)
xlim([0,2000])
title('fft spectral content')
xlabel('freq (Hz)')
ylabel('fa')
semilogx(freq,fa1,'color','r','linewidth',2)
hold on
semilogx(freq,fa2,'color','g','linewidth',2)
hold on
semilogx(freq,fa3,'color','k')

%%
% 1.3b/c

%1)
t = linspace(0,0.1,1000)
sin1 = sin(500*2*pi*t)
b=0
sin2_0 = sin(2*pi*500*(t-b))
b=0.0001
sin2_1 = sin(2*pi*500*(t-b))
b=0.0002
sin2_2 = sin(2*pi*500*(t-b))
b=0.0004
sin2_4 = sin(2*pi*500*(t-b))
b=0.0008
sin2_8 = sin(2*pi*500*(t-b))
b=0.001
sin2_10 = sin(2*pi*500*(t-b))

%2)
sin_sum_0 = sin1+sin2_0
sin_sum_1 = sin1+sin2_1
sin_sum_2 = sin1+sin2_2
sin_sum_4 = sin1+sin2_4
sin_sum_8 = sin1+sin2_8
sin_sum_10 = sin1+sin2_10
    
%3)
figure(4)
plot(t,sin_sum_0,'color','k','DisplayName','dt = 0s')
hold on
plot(t,sin_sum_1,'color','m','DisplayName','dt = 100 microsec')
hold on
plot(t,sin_sum_2,'color','c','DisplayName','dt = 200 microsec')
hold on
plot(t,sin_sum_4,'color','r','DisplayName','dt = 400 microsec')
hold on
plot(t,sin_sum_8,'color','b','DisplayName','dt = 800 microsec')
hold on
plot(t,sin_sum_10,'color','y','DisplayName','dt = 1000 microsec')
title('sine plot')
xlabel('t(s)')
ylabel('y')
xlim([0,0.01])
ylim([-2.5,2.5])
legend

%%
%4)
figure(5)
pwelch(sin_sum_0)
title('pwelch dt = 0s')

figure(6)
pwelch(sin_sum_1)
title('pwelch dt = 100 microsec')

figure(7)
pwelch(sin_sum_2)
title('pwelch dt = 200 microsec')

figure(8)
pwelch(sin_sum_4)
title('pwelch dt = 400 microsec')

figure(9)
pwelch(sin_sum_8)
title('pwelch dt = 800 microsec')

figure(10)
pwelch(sin_sum_10)
title('pwelch dt = 1000 microsec')

%interference occurs at 1000 microseconds. This is when the power/frequency
%is minimal at a very negative value. Also in the plot you can see the
%sines (yellow) cancel eachother out. 

%%
%1.3d)
dt = 0 
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(11)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 0s')
%I expect the spectrum to be very broad, because there are lots of
%frequencies that are summed. 

%1.3e)
dt = 0.0001
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(12)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 100 microsec')

dt = 0.0002
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(13)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 200 microsec')

dt = 0.0004
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(14)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 400 microsec')

dt = 0.0008
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(15)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 800 microsec')

dt = 0.001
sinfreq1 = 0
sinfreq2 = 0
for freq = 100:100:10000
    sinfreq1 = sinfreq1 + sin(2*pi*freq*t)
    sinfreq2 = sinfreq2 + sin(2*pi*freq*(t-dt))
    sinfreq_sum = sinfreq1+sinfreq2
end
figure(16)
pwelch(sinfreq_sum)
title('pwelch, broad spectrum, dt = 1000 microsec')

%The spectra differ in the sense that they are not noisy anymore, but there
%occurs a sin^2-like pattern with peaks at certain frequencies. Also there
%occur more peaks when the shift is bigger. 