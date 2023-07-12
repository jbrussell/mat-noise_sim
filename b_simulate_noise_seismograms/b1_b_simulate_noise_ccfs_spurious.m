% Simulate ambient noise within a homogeneous, isotropic, acoustic half-space.
% Following course notes from Goran's class and Wapenaar et al. (2005)
%
% This version includes a spurious source resulting in an asymmetric noise
% correlation function.
%
% jbrussell - 7/2023

clear; close all;

addpath('../functions/');

dx = 10; % km
x = [-1000:dx:1000]; % km
y = [-1000:dx:1000]; % km
dt = 1; % sec
Ndays = 10;
t = [0:dt:24*60*60 * Ndays]; % [sec] daily time axis

% Location of virtual source (A)
x_Asrc = 100; % km
y_Asrc = -100; % km

% Location of receiver (B)
x_Brec = -100; % km
y_Brec = 100; % km

% Location of "source" ring (S)
N_excite = 1; % number of times to randomly excite the wavelet
dtheta = 0.1;
theta_S = [0:dtheta:360-dtheta];
theta_S(end) = linspace(90-0,90+0,1);
r_S = 1000*ones(size(theta_S)); % [km] radious of source ring
r_S(end) = linspace(5000-0,5000+0,1);
x_S = r_S.*sind(theta_S);
y_S = r_S.*cosd(theta_S);
amp_S = ones(size(x_S)); % amplitude of sources
% Make one spurious noise source that is stronger than the rest
amp_S(end) = 8;

% Medium properties
vel = 3.5; % [km/s] velocity of medium

% Define source type
source_type = 'ricker'; % 'ricker' or 'microseism'
% RICKER  (impulsive source)
f_cent = 1/8; % 1/100; % [1/s] dominant frequency of Ricker wavelet
% MICROSEISM  (continuous source)
fmin = 1/10; % minimum frequency to sum over
fmax = 1/3; % maximum frequency to sum over

figure(98);
box on; hold on;
plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
text(x_Asrc+50,y_Asrc,'A','fontsize',13);
plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
text(x_Brec+50,y_Brec,'B','fontsize',13);
scatter(x_S,y_S,amp_S*20,amp_S,'sk','linewidth',2,'MarkerFaceColor','k');
text(x_S(1),y_S(1)+75,'S','fontsize',13);
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('X (km)');
ylabel('Y (km)');

%% Generate wavefield at ring of sources S

% Build frequency axis
Nt = length(t);
Fs = 1./dt;
f = Fs*(0:(Nt/2))/Nt;

% Generate wavefields at A and B due to random noise at all Si
Si_A = zeros(size(t));
Si_B = zeros(size(t));
for isrc = 1:length(x_S)
    R_Si_A = sqrt((x_Asrc-x_S(isrc)).^2 + (y_Asrc-y_S(isrc)).^2); % [km] distance from Si to A
    R_Si_B = sqrt((x_Brec-x_S(isrc)).^2 + (y_Brec-y_S(isrc)).^2); % [km] distance from Si to B

    % Loop through and generate sources, Si
    for ii = 1:N_excite
        if strcmp(source_type,'ricker')
            % Ricker wavelet
            % Generate random wavelet start times
            tshift = max(t) .* rand(1,1);
            Si_A = Si_A + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_A,vel,f_cent);
            Si_B = Si_B + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_B,vel,f_cent);
        elseif strcmp(source_type,'microseism')
            freq = f(f>=fmin & f<=fmax);
            phi_rand = (2*pi)*rand(1,length(freq)); % random phase between [0,2*pi]
            Si_A = Si_A + amp_S(isrc) .* microseism_source(t,R_Si_A,vel,freq,phi_rand);
            Si_B = Si_B + amp_S(isrc) .* microseism_source(t,R_Si_B,vel,freq,phi_rand);
        else
            error('Source type must be ''ricker'' or ''microseism''');
        end
    end
end

% Taper waveform
Si_A = cos_taper(Si_A);
Si_B = cos_taper(Si_B);

% Calculate power spectra
fftA = fft(Si_A);
% fftA = spectrumwhiten_smooth(fftA,0.001);
P_A = abs(fftA/Nt);
P_A = P_A(1:Nt/2+1);
P_A(2:end-1) = 2*P_A(2:end-1);
fftB = fft(Si_B);
% fftB = spectrumwhiten_smooth(fftB,0.001);
P_B = abs(fftB/Nt);
P_B = P_B(1:Nt/2+1);
P_B(2:end-1) = 2*P_B(2:end-1);

% Calculate cross-correlation
cohsum = fftB.*conj(fftA) ./ abs(fftB) ./ abs(fftA);
ccf = real(ifft(cohsum,Nt)); % inverse FFT to get time domain
ccf = fftshift(ccf); % rearrange values as [-lag lag]
ccf = detrend(ccf);
ccf = cos_taper(ccf);
time = ([0:Nt-1]-floor(Nt/2))*dt;  % build lagtime vector for plotting
time = [time(time<0), time(time>=0)];

% Expected time between A and B
R_A_B = sqrt((x_Asrc-x_Brec).^2 + (y_Asrc-y_Brec).^2); % [km] distance from A to B
t_A_B = R_A_B ./ vel;

% Theoretical prediction of ccf from Goran's class notes
ccf_pre = real(1 ./ sqrt(t_A_B.^2 - time.^2));

% Theoretical prediction of spurious source arrival time
azi_A_B = angle((y_Brec-y_Asrc) + 1i*(x_Brec-x_Asrc))*180/pi;
azi_S_A = angle((y_Asrc-y_S(end)) + 1i*(x_Asrc-x_S(end)))*180/pi;
t_A_B_spurious = t_A_B .* cosd(azi_A_B-azi_S_A);


%% Plot summary

figure(99); clf;
set(gcf,'position',[348         179        1040         799],'color','w');

subplot(3,2,1); box on; hold on;
plot(t,Si_A,'-g');
title('A');
set(gca,'fontsize',15,'linewidth',1.5)
xlabel('Time (s)');

subplot(3,2,3); box on; hold on;
plot(1./f,P_A,'-g');
plot(1./f,smooth(P_A,20),'-k','linewidth',2);
plot(1./f_cent*[1 1],[min(P_A) max(P_A)],'--r','linewidth',1.5);
% xlim([0 500]);
set(gca,'fontsize',15,'linewidth',1.5,'xscale','log','yscale','log')
xlabel('Period (s)');

subplot(3,2,2); box on; hold on;
plot(t,Si_B,'-b');
title('B');
set(gca,'fontsize',15,'linewidth',1.5)
xlabel('Time (s)');

subplot(3,2,4); box on; hold on;
plot(1./f,P_B,'-b');
plot(1./f,smooth(P_B,20),'-k','linewidth',2);
plot(1./f_cent*[1 1],[min(P_B) max(P_B)],'--r','linewidth',1.5);
% xlim([0 500]);
set(gca,'fontsize',15,'linewidth',1.5,'xscale','log','yscale','log')
xlabel('Period (s)');

subplot(3,2,[5 6]); box on; hold on;
plot(time,ccf ./ max(ccf),'-k','linewidth',2);
plot(time,ccf_pre ./ max(ccf_pre),'-r','linewidth',1);
plot(t_A_B_spurious*[1 1],[0 1],'--r','linewidth',1);
% plot(t_A_B*[1 1],[0 1],'--g','linewidth',1.5);
% plot(-t_A_B*[1 1],[0 1],'--g','linewidth',1.5);
set(gca,'fontsize',15,'linewidth',1.5)
xlim([-200 200]);
xlabel('Lag Time (s)');

save2pdf('noise_ccf_summary_spurious.pdf',99,300);
