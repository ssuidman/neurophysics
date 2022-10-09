%% Sweep signal
load('singlesweep.mat'); 
% load('HRTF.mat'); 
load('HRTF0119.mat');
t = linspace(0,length(sweep1)/fs,length(sweep1));
% We can measure until fs=10^5 Hz from this
sweep_fft = fft(sweep1,1024); % this is a good resolution for the FREQ (signal now happens to be the same), also can choose 2^n, but bigger n takes more computing power
f = linspace(1,fs,1024);
figure(1);
clf; % clean figure
subplot(3,1,1);
plot(t,sweep1); % waveform
xlabel('t(s)');
ylabel('amplitude');
title('sweep signal');
subplot(3,1,2);
semilogx(f,real(sweep_fft)); % magnitude
xlabel('f(Hz)');
ylabel('amplitude');
title('sweep spectrum magnitude');
subplot(3,1,3);
semilogx(f,imag(sweep_fft)); % phase
xlabel('f(Hz)');
ylabel('amplitude');
title('sweep spectrum phase');

% The sweep contains (almost) all frequencies and that is what you want if
% you measure the HRTF


%%
figure(2);
y = elSweep(:,el==0);
t = linspace(0,length(y)/fs,length(y));
plot(t,elSweep(:,el==0));
xlim([t(1) t(end)]);
xlabel('t(s)');
ylabel('amplitude');
title('waveform of central speaker');
shg;


%% Fourier of H

% central_sweep = elSweep(:,el==0)';
% t = linspace(0,length(central_sweep)/fs,length(central_sweep));
% sweep_fft = fft(central_sweep,1024); % don't know if this is right 

sweep_fft = fft(sweep1,1024);
f = linspace(0,fs/2,1024);

N = size(elSweep,2);
HRTF = zeros(N,1024);

figure(3);
clf;
for i = 1:N
    y = elSweep(:,i)';
    t = linspace(0,length(y)/fs,length(y));
    y_fft = fft(y,1024); % this is a good resolution for the FREQ (signal now happens to be the same), also can choose 2^n, but bigger n takes more computing power
    
    H = y_fft; 
%     H = y_fft./sweep_fft;
    HRTF(i,:) = abs(H);
    
    semilogx(f(1:512),log(abs(H(1:512))));
    hold on;
end
xlim([3000 12000]);
xticks([1000 2000 4000 8000 16000]);
xlabel('f(Hz)');
ylabel('amplitude(dB)');
title('HRTF');
shg;



%% 
DTF = getdtf(HRTF',Fs);
f_513 = f(1:513);

figure(5);
clf; % clean figure

subplot(2,1,1);
semilogx(f_513,DTF');
xlim([3000 12000]);
xlabel('f(Hz)');
ylabel('amplitude(db)');
title('DTF');

subplot(2,1,2);
imagesc(f_513,el,DTF');
colorbar; 
colormap jet; 
caxis([-15 10]);
set(gca,'YDir','normal');  
xlabel('f(Hz)');
ylabel('elevation (degrees)');
title('DTF');

