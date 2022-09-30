%%

%1.2a
t = linspace(-10*10^(-3),90*10^(-3),1000)

%1.2b
y1 = 1.5*sin(2*pi*150*t+pi/6)
%You expect to see about 150*2=130 peaks in one second (top and bottom of
%sine). In 100 ms = 0.1 second you expect to see 30 peaks (top and bottom).
figure(1)
plot(t,y1,'color','r','linewidth',2)

%1.2c
xlim([-0.01,0.09])
xlabel('t(s)')
ylabel('y')
title('plot of sinusoids')

%1.2d
hold on
y2 = 0.6*sin(2*pi*80*t+pi/3)
plot(t,y2,'color','g','linewidth',2)

%1.2e
hold on
y3 = y1+y2
plot(t,y3,'color','k')

%%
%1.2f
%As an example the two functions y1 and y2 can be plotted this way with plotsine(freq,phase,amplitude,colour):

figure(2)
xlim([-0.01,0.09])
ylim([-1.8,1.8])
plotsine(150,pi/6,1.5,'r')
hold on
plotsine(80,pi/3,0.5,'g')
title('plot out of plotsine function')
xlabel('t(s)')
ylabel('y')