%f

function sine = plotsine(freq,phase,amplitude,colour)
t = linspace(-10*10^(-3),90*10^(-3),1000)
sine = amplitude*sin(freq*2*pi*t+phase)

plot(t,sine,'color',colour,'linewidth',2)
xlabel('t')
ylabel('y')
title('plot of sinusoids')
end