% Simulate ambient noise within a homogeneous, isotropic, acoustic half-space 
% for random sources that form a finite-width donut and include a spurious
% source.
%
% jbrussell - 7/2023

clear; close all;

addpath('../functions/');

dx = 10; % km
x = [-1000:dx:1000]; % km
y = [-1000:dx:1000]; % km
dt = 1; % sec
Ndays = 10;
t = [0:dt:60*60]; % [sec] daily time axis

% Location of virtual source (A)
x_Asrc = 100; % km
y_Asrc = -100; % km

% Location of receiver (B)
x_Brec = -100; % km
y_Brec = 100; % km

% Location of "source" ring (S)
N_sources = 5000; % Number of sources
N_excite = 1; % number of times to randomly excite each Source each hour
dr_S = 200; % [km] width of annulus, S
r_S = 1000 + (dr_S*2*rand(N_sources,1)' - dr_S); %1000; % [km] radius of source ring
% dtheta = 0.1;
theta_S = 360*rand(N_sources,1)'; %[0:dtheta:360-dtheta];
% r_S(theta_S>90 & theta_S<180) = []; % remove sector of sources
% theta_S(theta_S>90 & theta_S<180) = []; % remove sector of sources
amp_S = ones(size(theta_S)); % amplitude of sources

% Set up sprious noise source as last element
r_S(end) = 5000; %linspace(5000-0,5000+0,1);
theta_S(end) = 90; %linspace(90-0,90+0,1);
amp_S(end) = 5*2;

x_S = r_S.*sind(theta_S);
y_S = r_S.*cosd(theta_S);

% Medium properties
vel = 3.5; % [km/s] velocity of medium

% Define source type
source_type = 'ricker'; % 'ricker' | 'microseism' | 'lossy_membrane'
% RICKER (impulsive source)
f_cent = 1/8; % 1/100; % [1/s] dominant frequency of Ricker wavelet
% MICROSEISM (continuous source)
fmin = 1/10; % minimum frequency to sum over
fmax = 1/3; % maximum frequency to sum over
% LOSSY MEMBRANE (attenuating membrane)
fmin;
fmax;
alpha = 1e-3; % [1/km] attenuation coefficient


figure(98); clf;
box on; hold on;
plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
text(x_Asrc+50,y_Asrc,'A','fontsize',13);
plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
text(x_Brec+50,y_Brec,'B','fontsize',13);
scatter(x_S,y_S,amp_S*5,amp_S,'o','filled');
text(x_S(1),y_S(1)+75,'S','fontsize',13);
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('X (km)');
ylabel('Y (km)');
cb = colorbar;
ylabel(cb,'Source Amplitude');
colormap(viridis);

%% Generate wavefield at ring of sources S

% Build frequency and time axes
Nt = length(t);
Fs = 1./dt;
f = Fs*(0:(Nt/2))/Nt;
time = ([0:Nt-1]-floor(Nt/2))*dt;  % build lagtime vector for plotting
time = [time(time<0), time(time>=0)];

% Expected time between A and B
R_A_B = sqrt((x_Asrc-x_Brec).^2 + (y_Asrc-y_Brec).^2); % [km] distance from A to B
t_A_B = R_A_B ./ vel;

% Theoretical prediction of ccf from Goran's class notes
ccf_pre = real(1 ./ sqrt(t_A_B.^2 - time.^2));
A = 1;
J0_pre = besselj(0,2*pi*f*t_A_B) * A;

% Theoretical prediction of spurious source arrival time
azi_A_B = angle((y_Brec-y_Asrc) + 1i*(x_Brec-x_Asrc))*180/pi;
azi_S_A = angle((y_Asrc-y_S(end)) + 1i*(x_Asrc-x_S(end)))*180/pi;
t_A_B_spurious = t_A_B .* cosd(azi_A_B-azi_S_A);

cohsum = zeros(size(t));
ccf_auto = zeros(size(t));
ihrs_total = 0;
ccf_all = {};
cohsum_all = {};
ccf_misfit = [];
for iday = 1:Ndays
    for ihr = 1 : (24/(max(t)/60/60))
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
                    Si_A = Si_A + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_A(isrc),vel,f_cent);
                    Si_B = Si_B + amp_S(isrc) .* ricker_wavelet(t-tshift,R_Si_B(isrc),vel,f_cent);
                elseif strcmp(source_type,'microseism')
                    % Continuous microseism source
                    freq = f(f>=fmin & f<=fmax);
                    phi_rand = (2*pi)*rand(1,length(freq)); % random phase between [0,2*pi]
                    Si_A = Si_A + amp_S(isrc) .* microseism_source(t,R_Si_A(isrc),vel,freq,phi_rand);
                    Si_B = Si_B + amp_S(isrc) .* microseism_source(t,R_Si_B(isrc),vel,freq,phi_rand);
                elseif strcmp(source_type,'lossy_membrane')
                    % 2-D Lossy Membrane response of Magrini & Boschi (2021)
                    phi_rand = (2*pi)*rand(1,1); % random phase between [0,2*pi]
                    Si_A = Si_A + amp_S(isrc) .* lossy_membrane(t,f,R_Si_A(isrc),vel,alpha,phi_rand);
                    Si_B = Si_B + amp_S(isrc) .* lossy_membrane(t,f,R_Si_B(isrc),vel,alpha,phi_rand);
                else
                    error('Source type must be ''ricker'' | ''microseism'' | ''lossy_membrane''');
                end
            end
        end
        
        % Taper waveform
        Si_A = cos_taper(Si_A);
        Si_B = cos_taper(Si_B);

        % Calculate power spectra
        fftA = fft(Si_A);
        % fftA = spectrumwhiten_smooth(fftA,0.001);
%         P_A = abs(fftA/Nt);
%         P_A = P_A(1:Nt/2+1);
%         P_A(2:end-1) = 2*P_A(2:end-1);
        fftB = fft(Si_B);
        % fftB = spectrumwhiten_smooth(fftB,0.001);
%         P_B = abs(fftB/Nt);
%         P_B = P_B(1:Nt/2+1);
%         P_B(2:end-1) = 2*P_B(2:end-1);

        % Calculate cross-correlation (A --> B Causal; B --> A Acausal)
        cohsum = cohsum + fftB.*conj(fftA) ./ abs(fftB) ./ abs(fftA);
        ccf = real(ifft(cohsum,Nt)); % inverse FFT to get time domain
        ccf = fftshift(ccf); % rearrange values as [-lag lag]
        ccf = detrend(ccf);
        ccf = cos_taper(ccf);
        
        % Auto-correlation
        ccf_auto = ccf_auto + xcorr(Si_A,Si_A,floor(length(t)/2));
        
        ihrs_total = ihrs_total + 1;
        ccf_all{ihrs_total} = ccf;
        cohsum_all{ihrs_total} = cohsum;
        ccf_misfit(ihrs_total) = sum((ccf_all{ihrs_total} ./ max(ccf_all{ihrs_total}) - ccf_pre ./ max(ccf_pre)).^2) / length(ccf_pre);
        
        if 1 && mod(ihrs_total,12)==0
            figure(99); clf;
            set(gcf,'position',[2         206        1730         814]);
            
            subplot(2,2,1); box on; hold on;
            plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
            text(x_Asrc+50,y_Asrc,'A','fontsize',13);
            plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
            text(x_Brec+50,y_Brec,'B','fontsize',13);
            scatter(x_S,y_S,amp_S*5,amp_S,'o','filled');
            text(x_S(1),y_S(1)+75,'S','fontsize',13);
            axis square; axis equal;
            set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
            xlabel('X (km)');
            ylabel('Y (km)');
            cb = colorbar;
            ylabel(cb,'Source Amplitude');
            colormap(viridis);
            
            subplot(2,2,3); box on; hold on;
            plot([1:ihrs_total],ccf_misfit,'-b','linewidth',2);
            plot(ihrs_total,ccf_misfit(ihrs_total),'ok','linewidth',1,'markersize',12,'MarkerFaceColor','b');
            set(gca,'fontsize',15,'linewidth',1.5);
            xlabel('Data length (hr)');
            ylabel('Misfit');
            
%             % Filter CCF before plotting
%             coperiod = [5 20];
%             costap_wid = 0.2; % 0 => box filter; 1 => Hann window
%             [ ccf_filtered ] = tukey_filt( fft(fftshift(ccf_all{ihrs_total})),coperiod,dt,costap_wid );
%             ccf_all{ihrs_total} = fftshift(real(ifft(ccf_filtered)));

            subplot(2,2,2); box on; hold on;
%             ccf_pre_conv = conv(ccf_pre,ccf_auto,'same');
            plot(time,ccf_all{ihrs_total} ./ max(ccf_all{ihrs_total}),'-k','linewidth',2);
            plot(time,ccf_pre ./ max(ccf_pre),'-r','linewidth',1);
%             plot(time,ccf_pre_conv ./ max(ccf_pre_conv),'--g','linewidth',1);
            plot(t_A_B_spurious*[1 1],[0 1],'--r','linewidth',1);
            set(gca,'fontsize',15,'linewidth',1.5)
            xlim([-200 200]);
            title([num2str(ihrs_total),' hrs']);
            xlabel('Lag Time (s)');
            ylabel('CCF_{A,B}');

            P_ccf = abs(cohsum_all{ihrs_total}/Nt);
            P_ccf = P_ccf(1:Nt/2+1);
            P_ccf(2:end-1) = 2*P_ccf(2:end-1);

%             subplot(2,1,2); box on; hold on;
%             plot(1./f,P_ccf,'-b');
%             plot(1./f,smooth(P_ccf,20),'-k','linewidth',2);
%             plot(1./f_cent*[1 1],[min(P_ccf) max(P_ccf)],'--r','linewidth',1.5);
%             % xlim([0 500]);
%             set(gca,'fontsize',15,'linewidth',1.5,'xscale','log','yscale','log')
%             xlabel('Period (s)');

            subplot(2,2,4); box on; hold on;
            J0 = real(cohsum_all{ihrs_total}(1:Nt/2+1));
            plot(f,J0 ./ rms(J0),'-b');
            plot(f,smooth(J0,10) ./ rms(J0),'-k','linewidth',2);
            plot(f,J0_pre ./ rms(J0_pre),'-r','linewidth',1);
            plot(f_cent*[1 1],[0 1],'--r','linewidth',1.5);
            % xlim([0 500]);
            set(gca,'fontsize',15,'linewidth',1.5,'xscale','linear','yscale','linear')
            xlabel('Freq (Hz)');
            ylabel('J_0');
            xlim([0 max(f)]);

            drawnow;
        end
    end
end

%%
clear OUTmovie

is_movieout = 1;

if is_movieout
        
    OUTmovie(length(1:1:ihrs_total)) = struct('cdata',[],'colormap',[]);
    OUTit = 0;
    misfit = [];
    for ii = 1:1:ihrs_total
    
        figure(3); clf;
        set(gcf,'position',[2         206        1730         814],'color','w');
        
        subplot(2,2,1); box on; hold on;
        plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
        text(x_Asrc+50,y_Asrc,'A','fontsize',13);
        plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
        text(x_Brec+50,y_Brec,'B','fontsize',13);
        scatter(x_S,y_S,amp_S*5,amp_S,'o','filled');
        text(x_S(1),y_S(1)+75,'S','fontsize',13);
        axis square; axis equal;
        set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
        xlabel('X (km)');
        ylabel('Y (km)');
        cb = colorbar;
        ylabel(cb,'Source Amplitude');
        colormap(viridis);
        
        subplot(2,2,3); box on; hold on;
        plot([1:ii],ccf_misfit(1:ii),'-b','linewidth',2);
        plot(ii,ccf_misfit(ii),'ok','linewidth',1,'markersize',12,'MarkerFaceColor','b');
        set(gca,'fontsize',15,'linewidth',1.5);
        xlim([0 length(ccf_misfit)]);
        ylim([0 max(ccf_misfit)]);
        xlabel('Data length (hr)');
        ylabel('Misfit');
        
        subplot(2,2,2); box on; hold on;
%         ccf_pre_conv = conv(ccf_pre,ccf_auto,'same');
        plot(time,ccf_all{ii} ./ max(ccf_all{ii}),'-k','linewidth',2);
        plot(time,ccf_pre ./ max(ccf_pre),'-r','linewidth',1);
%         plot(time,ccf_pre_conv ./ max(ccf_pre_conv),'--g','linewidth',1);
        plot(t_A_B_spurious*[1 1],[0 1],'--r','linewidth',1);
        set(gca,'fontsize',15,'linewidth',1.5)
        xlim([-200 200]);
        ylim([-0.3 1.1]);
        title([num2str(ii),' hrs']);
        xlabel('Lag Time (s)');
        ylabel('CCF_{A,B}');

        P_ccf = abs(cohsum_all{ii}/Nt);
        P_ccf = P_ccf(1:Nt/2+1);
        P_ccf(2:end-1) = 2*P_ccf(2:end-1);

%         subplot(2,1,2); box on; hold on;
%         plot(1./f,P_ccf,'-b');
%         plot(1./f,smooth(P_ccf,20),'-k','linewidth',2);
%         plot(1./f_cent*[1 1],[min(P_ccf) max(P_ccf)],'--r','linewidth',1.5);
%         % xlim([0 500]);
%         set(gca,'fontsize',15,'linewidth',1.5,'xscale','log','yscale','log')
%         xlabel('Period (s)');
        
        subplot(2,2,4); box on; hold on;
        J0 = real(cohsum_all{ii}(1:Nt/2+1));
        plot(f,J0 ./ rms(J0),'-b');
        plot(f,smooth(J0,10) ./ rms(J0),'-k','linewidth',2);
        plot(f,J0_pre ./ rms(J0_pre),'-r','linewidth',1);
        plot(f_cent*[1 1],[0 1],'--r','linewidth',1.5);
        % xlim([0 500]);
        set(gca,'fontsize',15,'linewidth',1.5,'xscale','linear','yscale','linear')
        xlabel('Freq (Hz)');
        ylabel('J_0');
        xlim([0 max(f)]);

        drawnow;

        OUTit = OUTit + 1;
        rect = get(gcf,'Position'); rect(1:2) = [0 0]; 
        OUTmovie(OUTit) = getframe(gcf,rect);
    end
    
%     vid = VideoWriter(['noise_ccf_build_Sdonut_spurious','.avi']);
    vid = VideoWriter(['noise_ccf_build_Sdonut_spurious'],'MPEG-4');
    vid.FrameRate = 10;
    open(vid);
    OUTmovie(1) = [];
    writeVideo(vid,OUTmovie);
    close(vid);
end
