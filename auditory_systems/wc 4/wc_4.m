% load('test_data.mat', 'Data');
load('NH-0117-22-10-14-data.mat', 'Data');

%% test figure
stimaz = Data(:,1); % stimulus azimuth in first column
resaz = Data(:,5); % response azimuth in 5th column

figure(1);
plot(stimaz,resaz,'ko'); % plot responses as a function of azimuth
xlabel('stimulus azimuth');
ylabel('response azimuth');
xlim([-90 90]);
ylim([-90 90]);
hold on;
plot([-90 90],[-90 90],':');

%% (2,3)-plot for elevation and frequency bands

freq = Data(:,3);

% azimuth BB
stimaz_BB = Data(freq==1, 1); % stimulus azimuth in first column, where freq==1 means broadband sounds
resaz_BB = Data(freq==1, 5); % response azimuth in 5th column
subplot(2,3,1);
plot(stimaz_BB,resaz_BB,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('azimuth, BB: 0.2–20 kHz');
hold on;
plot([-90 90],[-90 90]);

% azimuth HP
stimaz_HP = Data(freq==2, 1); % stimulus azimuth in first column, where freq==1 means broadband sounds
resaz_HP = Data(freq==2, 5); % response azimuth in 5th column
subplot(2,3,2);
plot(stimaz_HP,resaz_HP,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('azimuth, HP: 3-20 kHz');
hold on;
plot([-90 90],[-90 90]);

% azimuth LP
stimaz_LP = Data(freq==3, 1); % stimulus azimuth in first column, where freq==1 means broadband sounds
resaz_LP = Data(freq==3, 5); % response azimuth in 5th column
subplot(2,3,3);
plot(stimaz_LP,resaz_LP,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('azimuth, LP: 0.2-1.5 kHz');
hold on;
plot([-90 90],[-90 90]);

% elevation BB
stimel_BB = Data(freq==1, 2); % stimulus azimuth in first column, where freq==1 means broadband sounds
resel_BB = Data(freq==1, 6); % response azimuth in 5th column
subplot(2,3,4);
plot(stimel_BB,resel_BB,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('elevation, BB: 0.2–20 kHz');
hold on;
plot([-90 90],[-90 90]);

% elevation HP
stimel_HP = Data(freq==2, 2); % stimulus azimuth in first column, where freq==1 means broadband sounds
resel_HP = Data(freq==2, 6); % response azimuth in 5th column
subplot(2,3,5);
plot(stimel_HP,resel_HP,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('elevation, HP: 3.0–20 kHz');
hold on;
plot([-90 90],[-90 90]);

% elevation LP
stimel_LP = Data(freq==3, 2); % stimulus azimuth in first column, where freq==1 means broadband sounds
resel_LP = Data(freq==3, 6); % response azimuth in 5th column
subplot(2,3,6);
plot(stimel_LP,resel_LP,'ko'); % plot responses as a function of azimuth
xlim([-90 90]);
ylim([-90 90]);
xlabel('stimulus');
ylabel('response');
title('elevation, LP: 0.2-1.5 kHz');
hold on;
plot([-90 90],[-90 90]);



%% slope and offset for 6 possibilities

B_az_BB = regstats(resaz_BB,stimaz_BB,'linear');
slope_az_BB = B_az_BB.beta(2);
offset_az_BB = B_az_BB.beta(1);

B_az_HP = regstats(resaz_HP,stimaz_HP,'linear');
slope_az_HP = B_az_HP.beta(2);
offset_az_HP = B_az_HP.beta(1);

B_az_LP = regstats(resaz_LP,stimaz_LP,'linear');
slope_az_LP = B_az_LP.beta(2);
offset_az_LP = B_az_LP.beta(1);

B_el_BB = regstats(resel_BB,stimel_BB,'linear');
slope_el_BB = B_el_BB.beta(2);
offset_el_BB = B_el_BB.beta(1);

B_el_HP = regstats(resel_HP,stimel_HP,'linear');
slope_el_HP = B_el_HP.beta(2);
offset_el_HP = B_el_HP.beta(1);

B_el_LP = regstats(resel_LP,stimel_LP,'linear');
slope_el_LP = B_el_LP.beta(2);
offset_el_LP = B_el_LP.beta(1);

display(slope_az_BB);
display(slope_az_HP);
display(slope_az_LP);

display(slope_el_BB);
display(slope_el_HP);
display(slope_el_LP);

display(offset_az_BB);
display(offset_az_HP);
display(offset_az_LP);

display(offset_el_BB);
display(offset_el_HP);
display(offset_el_LP);







