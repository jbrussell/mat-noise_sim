% Simulate a perfect ring of sources
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
r_S = 1000; % [km] radious of source ring
N_excite = 1; % number of times to randomly excite each Source each hour
dtheta = 0.1;
theta_S = [0:dtheta:360-dtheta];
x_S = r_S*sind(theta_S);
y_S = r_S*cosd(theta_S);
amp_S = ones(size(x_S)); % amplitude of sources

% Gaussian/Ricker wavelet properties
% c = 3.5; % [km/s] phase velocity
grv = 3.5; % [km/s] group velocity
% sig = 40; % km;
f_cent = 1/8; % 1/100; % [1/s] freq


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
t_A_B = R_A_B ./ grv;

% Theoretical prediction of ccf from Goran's class notes
ccf_pre = real(1 ./ sqrt(t_A_B.^2 - time.^2));
A = 1;
J0_pre = besselj(0,2*pi*f*t_A_B) * A;

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

            % Ricker wavelet
            % Generate random wavelet start times
            tshift = sort(0 + (max(t)-0) .* rand(N_excite,1));
            % Loop through and generate Ricker wavelets for Si
            for ii = 1:N_excite
                Si_A = Si_A + amp_S(isrc) .* ricker_wavelet(t-tshift(ii),R_Si_A,grv,f_cent);
                Si_B = Si_B + amp_S(isrc) .* ricker_wavelet(t-tshift(ii),R_Si_B,grv,f_cent);
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
%         ccf_auto = ccf_auto + xcorr(Si_A,Si_A,floor(length(t)/2));
        
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
            axis square;
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

            subplot(2,2,2); box on; hold on;
%             ccf_pre_conv = conv(ccf_pre,ccf_auto,'same');
            plot(time,ccf_all{ihrs_total} ./ max(ccf_all{ihrs_total}),'-k','linewidth',2);
            plot(time,ccf_pre ./ max(ccf_pre),'-r','linewidth',1);
%             plot(time,ccf_pre_conv ./ max(ccf_pre_conv),'--g','linewidth',1);
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
    
    vid = VideoWriter(['noise_ccf_build_Sring','.avi']);
    vid.FrameRate = 10;
    open(vid);
    OUTmovie(1) = [];
    writeVideo(vid,OUTmovie);
    close(vid);
end
