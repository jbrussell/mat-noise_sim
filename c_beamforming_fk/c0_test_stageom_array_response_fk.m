% Calculate Array Response Function (ARF). The ARF depends only on the 
% geometry of the seismic array and not on the data.
%
% See Eq. 2.16 of "Seismic Ambient Noise" by Nakata, Gualtieri, & Fichtner
%
% jbrussell - 2/2024

addpath('../functions/');

clear; close all;
rng default % for reproducibility

%========== Define plane wave properties for Array Response Function ==========%

f_p = 1/8; % [Hz] frequency of plane wave
az_p = 0; %[deg] propagation azimuth. 0=North propagation (i.e., arriving from south)
slow_p = 0; %1/3.5; %[s/km] slowness of plane wave (zero gives vertically incident wave)
% Slowness vector pointing from source N toward array
%                  x           y
s_p = slow_p*[sind(az_p); cosd(az_p)]; 

%======================= BUILD ARRAY =======================%

% Origin in lat lon
olat = 35;
olon = -105;
deg = 0.5; % station spacing

% Define synthetic array of stations
[Y_stas, X_stas] = meshgrid(olat+[-1:deg:1], olon+[-1:deg:1]);
X_stas = X_stas(:);
Y_stas = Y_stas(:);
% Perturb station locations slightly
d_deg = km2deg(15);
X_stas = X_stas + d_deg*2*rand(length(X_stas),1)-d_deg;
Y_stas = Y_stas + d_deg*2*rand(length(Y_stas),1)-d_deg;

%======================= BEAMFORMING PARAMETERS =======================%

% Data weights and SNR threshold
is_azimuthal_weight = 1; % apply azimuthal weights to downweight redundant azimuths?
dazi = 10; % [deg] bin width for azimuthal homogenization
is_dist_weight = 1; % down weight shorter station separations

% Periods to average over
per_min = 1/f_p; %5; % [sec] minimum period
per_max = 1/f_p; %10; % [sec] maximum period
Npers = 1;  %30; % number of periods to consider
f_vec = linspace(1/per_max,1/per_min,Npers);
per_vec = 1./f_vec;

% Slowness values to search over
s_min = 0; %1/5; % s/km
s_max = 1/2.5; % s/km
Nslow = 100;
s_vec = linspace(s_min,s_max,Nslow);

% Back-azimuth values to search over
Nbaz = 360;
baz_vec = linspace(0,360,Nbaz);

% =======================================================================

% Plot station geometry

figure(98); clf;
box on; hold on;
load coastlines
plot(coastlon,coastlat,'-b');
plot(X_stas,Y_stas,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5)
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('Lon (X)');
ylabel('Lat (Y)');
xlim([min(X_stas)-1 max(X_stas)+1]);
ylim([min(Y_stas)-1 max(Y_stas)+1]);

%% Determine data weights based on azimuth (downweight common azimuths)

nsta = length(X_stas); % number of target stations to calculate for

% Loop through all stations and gather info
stainfo = [];
ii = 0;
for ista1 = 1:nsta
    sta1 = ['STA',num2str(ista1,'%.2d')];
    stainfo.stas{ista1} = sta1;
    for ista2 = ista1:nsta
        
        sta2 = ['STA',num2str(ista2,'%.2d')];
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end

        disp([sta1,'-',sta2]);
        
        r1 = distance(Y_stas(ista1),X_stas(ista1),Y_stas(ista2),X_stas(ista2),referenceEllipsoid('GRS80'))/1000;
        [~,az] = distance(Y_stas(ista1),X_stas(ista1),Y_stas(ista2),X_stas(ista2),referenceEllipsoid('GRS80'));
                
        ii = ii + 1;
        stainfo.r(ii) = r1;
        stainfo.az(ii) = az;
        stainfo.sta1{ii} = sta1;
        stainfo.sta2{ii} = sta2;
    end
end

stainfo.w = ones(size(stainfo.r));

% Down weight redundant azimuths
if is_azimuthal_weight
    az_bins = [0:10:360];
    az_weights = ones(size(az_bins));
    for ibin = 1:length(az_bins)-1
        I = find(stainfo.az>=az_bins(ibin) & stainfo.az<az_bins(ibin+1));
        az_weights(ibin) = 1/length(I);
        stainfo.w(I) = stainfo.w(I) * (1/length(I));
    end
end

% Down weight shorter station separations
if is_dist_weight
    dist_weights = (stainfo.r/max(stainfo.r));
    stainfo.w = stainfo.w .* (stainfo.r/max(stainfo.r));
%     stainfo.w(stainfo.r<200) = 0;
end

%% Do beamforming

%%% --- Loop through station 1 --- %%%
Pf = zeros(Nbaz,Nslow,Npers);
for ista1= 1:nsta   
    sta1=stainfo.stas{ista1};
    
    %%% --- Loop through station 2 --- %%%
    for ista2 = ista1: nsta % length(v_sta)
        sta2 = stainfo.stas{ista2};
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
        
        % Get weight for station pair
        Ipair = find(strcmp(stainfo.sta1,sta1) & strcmp(stainfo.sta2,sta2));
        w = stainfo.w(Ipair);
        if isempty(w)
            continue
        end
        
        disp([sta1,'-',sta2]);
        
        r1 = stainfo.r(Ipair);
        az = stainfo.az(Ipair);
        
        % Displacement vector pointing from station 1 to station 2
        %          x             y
        r_v = r1*[sind(az); cosd(az)];
        
        for iper = 1:Npers
            per = per_vec(iper);
            omega = 2*pi ./ per;
            
            for islow = 1:Nslow
                slow = s_vec(islow);
                
                for ibaz = 1:Nbaz
                    baz = baz_vec(ibaz);
                    
                    % Slowness vector pointing from source N toward array
                    %                  x             y
                    s_v = slow*[sind(baz-180); cosd(baz-180)]; 
                    
                    % slowness direction dotted onto interstation path
                    sdotr = (s_v-s_p)'*r_v;
                    
                    % Estimate power
                    Pf(ibaz,islow,iper) = Pf(ibaz,islow,iper) + w * exp(1i*omega * sdotr);
                end
            end
        end
%         figure(1); clf;
%         % scatter(x(:),y(:),10,P_abs(:));
%         % contourf(x,y,P_abs,'LineColor','none'); axis equal;
%         % polarscatter(baz_mat(:)*pi/180,s_mat(:),10,P_abs(:));
%         [h,c] = polarPcolor(s_vec,baz_vec,abs(P),'Nspokes',9);
%         colormap(viridis);
    end
end

% Sum over frequencies
P = sum(Pf,3);

%% Plot beam

P_abs = abs(P);
P_abs = P_abs / max(P_abs(:));
P_abs = 10*log10(P_abs);

figure(1); clf;

subplot(1,2,1);
box on; hold on;
load coastlines
plot(coastlon,coastlat,'-b');
plot(X_stas,Y_stas,'ok','markerfacecolor','r','markersize',8,'linewidth',1.5)
axis square; axis equal;
set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
xlabel('Lon (X)');
ylabel('Lat (Y)');
xlim([min(X_stas)-1 max(X_stas)+1]);
ylim([min(Y_stas)-1 max(Y_stas)+1]);
title('Station Geometry');

subplot(1,2,2)
set(gcf,'position',[183         328        1254         527],'color','w');
[h,c] = polarPcolor(s_vec,baz_vec,P_abs,'Nspokes',9,'fontsize',13);
colormap(viridis);
c.LineWidth = 1.5;
ylabel(c,'Relative Power (dB)','fontsize',15);
set(gca,'fontsize',15,'linewidth',1.5)
caxis([prctile(P_abs(:),80) 0]);
titl = title([num2str(1/f_p),'s']);
titl.Position(2) = titl.Position(2) + 0.25;
hp = polar((az_p+180+90)*pi/180,slow_p/(s_max-s_min),'-or');
hp.LineWidth = 1.5;
hp.MarkerEdgeColor = 'w';
hp.MarkerFaceColor = 'r';

% save2pdf([figpath,'fk_array_response_',num2str(1/f_p),'s.pdf'],1,250);


