% Following Wapenaar et al. (2005) time-reversal approach
%
% The purpose is to give physical intuition behind the acausal (B --> A) 
% and causal (A --> B) arrivals in ambient noise correlation functions.  
% 
% This naive version generates a source field at A (causal), time reverses it 
% (acausal), and concatinates it to the front of the causal field.
%
% jbrussell - 7/2023

clear; close all;

addpath('../functions/');

dx = 10; % km
x = [-1000:dx:1000]; % km
y = [-1000:dx:1000]; % km
dt = 1; % sec
% t = [-100:dt:100]; % sec
t_causal = [0:dt:500];
[X_causal,Y_causal,T_causal] = meshgrid(x,y,t_causal); 

% Location of virtual source (A)
x_Asrc = 200; % km
y_Asrc = -500; % km

% Location of receiver (B)
x_Brec = -200; % km
y_Brec = 200; % km

% Location of "source" ring (S)
r_S = 750; % [km] radius of source ring
dtheta = 5; % [deg] spacing between sources
theta_S = [0:dtheta:360-dtheta];
x_S = r_S*sind(theta_S);
y_S = r_S*cosd(theta_S);
amp_S = ones(size(x_S)); % amplitude of sources

% Gaussian wavelet properties
% c = 3.5; % [km/s] phase velocity
grv = 3.5; % [km/s] group velocity
% sig = 40; % km;
f_cent = 1/50; %1/100; % [1/s] freq

%% Calculate virtual source field (A)

% Time slice for plotting
i_tslice = 750; % [sec]

% Distance from source
R_A = sqrt((X_causal-x_Asrc).^2 + (Y_causal-y_Asrc).^2);

% % construct gaussian wavelet wavefield
% omega = 2*pi*f_cent; % angular freq.
% k = omega ./ c; % wavenumber
% % Gaussian wavelet = sinusoid * gaussian
% A_causal = exp(1i*(k.*R_A - omega.*T_causal)) .* 1./(sig*sqrt(2*pi)).*exp(-(R_A-grv.*T_causal).^2 ./ (2*sig^2));
% A_causal = real(A_causal);

% Ricker wavelet
A_causal = ricker_wavelet(T_causal,R_A,grv,f_cent);

% Generate acasual field
A_acausal = flip(A_causal,3);
t_acausal = -1*flip(t_causal);

% Build acausal + causal field
A = cat(3,A_acausal,A_causal(:,:,2:end));
t = [t_acausal, t_causal(2:end)];
[X,Y,T] = meshgrid(x,y,t); 


figure(1); clf;

box on; hold on;
Xslice = X(:,:,i_tslice);
Yslice = Y(:,:,i_tslice);
Aslice = A(:,:,i_tslice);
scatter(Xslice(:),Yslice(:),10,Aslice(:),'filled');
plot(x_Asrc,y_Asrc,'oc','linewidth',2,'MarkerFaceColor','c');
text(x_Asrc+50,y_Asrc,'A');
plot(x_Brec,y_Brec,'og','linewidth',2,'MarkerFaceColor','g');
text(x_Brec+50,y_Brec,'B');
scatter(x_S,y_S,amp_S*20,amp_S,'sm','linewidth',2,'MarkerFaceColor','m');
colormap(redblue)
cb = colorbar;
caxis([-max(abs(cb.Limits)) max(abs(cb.Limits))]);

%% Evaluate wavefield at S and time reverse back to A

u_B = zeros(1,length(t));
u_A = zeros(1,length(t));
% A field "recorded" at S ring
for it = 1:length(t)
    u_B(it) = interp2(X(:,:,it),Y(:,:,it),A(:,:,it),x_Brec,y_Brec);
    u_A(it) = interp2(X(:,:,it),Y(:,:,it),A(:,:,it),x_Asrc,y_Asrc);
end

figure(2); clf;

subplot(2,1,1);
box on; hold on;
plot(t,u_A,'-c','linewidth',2);
title('u_A')
xlabel('Time (s)');

subplot(2,1,2);
box on; hold on;
plot(t,u_B,'-g','linewidth',2);
title('u_B');
xlabel('Time (s)');

%% Write movie

is_movieout = 1;

if is_movieout
    
    it_plot = [1:25:length(t)];
    
    OUTmovie(length(it_plot)) = struct('cdata',[],'colormap',[]);
    OUTit = 0;
    for i_tslice = it_plot
    
        figure(3); clf;
        set(gcf,'position',[205   220   971   798]);
        
        subplot(2,2,[1 3]);
        box on; hold on;
        Xslice = X(:,:,i_tslice);
        Yslice = Y(:,:,i_tslice);
        Aslice = A(:,:,i_tslice);
        scatter(Xslice(:),Yslice(:),10,Aslice(:),'filled');
        plot(x_Asrc,y_Asrc,'og','linewidth',2,'MarkerFaceColor','g');
        text(x_Asrc+50,y_Asrc,'A','fontsize',13);
        plot(x_Brec,y_Brec,'ob','linewidth',2,'MarkerFaceColor','b');
        text(x_Brec+50,y_Brec,'B','fontsize',13);
        scatter(x_S,y_S,amp_S*20,amp_S,'sk','linewidth',2,'MarkerFaceColor','k');
        text(x_S(1),y_S(1)+75,'S','fontsize',13);
        colormap(redblue)
        cb = colorbar;
        caxis([-max(abs(u_B))*0.9 max(abs(u_B))*0.9]);
        axis square;
        set(gca,'fontsize',15,'linewidth',1.5,'layer','top');
        title('$u(x,t) = G(x,x_A,t) + G(x,x_A,-t)$','Interpreter','latex');
        
        subplot(2,2,2);
        box on; hold on;
        plot(t(1:i_tslice),u_A(1:i_tslice),'-g','linewidth',2);
        plot(t(i_tslice),u_A(i_tslice),'ok','markerfacecolor','g','markersize',12,'linewidth',1);
        title('$u(x_A,t) = G(x_A,x_A,t) + G(x_A,x_A,-t)$','Interpreter','latex')
        xlabel('Time (s)');
        xlim([min(t) max(t)]);
        ylim([min(u_A) max(u_A)]);
        set(gca,'fontsize',15,'linewidth',1.5);

        subplot(2,2,4);
        box on; hold on;
        plot(t(1:i_tslice),u_B(1:i_tslice),'-b','linewidth',2);
        plot(t(i_tslice),u_B(i_tslice),'ok','markerfacecolor','b','markersize',12,'linewidth',1);
        title('$u(x_B,t) = G(x_B,x_A,t) + G(x_B,x_A,-t)$','Interpreter','latex')
        xlabel('Time (s)');
        xlim([min(t) max(t)]);
        ylim([min(u_B) max(u_B)]);
        set(gca,'fontsize',15,'linewidth',1.5);
        
        drawnow;

        OUTit = OUTit + 1;
        rect = get(gcf,'Position'); rect(1:2) = [0 0]; 
        OUTmovie(OUTit) = getframe(gcf,rect);
    end
    
    vid = VideoWriter(['time_reversal_naive','.avi']);
    vid.FrameRate = 10;
    open(vid);
    OUTmovie(1) = [];
    writeVideo(vid,OUTmovie);
    close(vid);
end
