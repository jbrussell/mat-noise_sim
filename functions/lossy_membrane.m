function [wavefield] = lossy_membrane(t,freq,r,vel,alpha,phi_rand)
% Calculate frequeny-domain Green's function and seismic response for a 2-D 
% lossy membrane following Magrini & Boschi (2021) GJI "Surface-Wave 
% Attenuation From Seismic Ambient Noise: Numerical Validation and 
% Application"
%
% t: time vector for plotting [s]
% faxis: frequency vector [Hz]
% r: distance from source [km]
% vel: velocity of medium [km/s]
% alpha: attenuation coefficient [1/km]
% phi_rand: random uniform phase between [0,2*pi] rad
%
% jbrussell - 7/2023

% Angular frequency
omega = 2*pi*freq;

% % Approximate (alpha << omega/c) frequency-domain Green's function (eq 2 of Magrini & Boschi, 2021)
% prefac = -1i / (4*sqrt(2*pi)*vel^2);
% arg = (omega * r) / vel;
% % Zero-order Hankel function of the 2nd kind
% gf_f = prefac * besselh(0,2,arg) * exp(-alpha * r);

% Exact frequency-domain Green's function (eq 1 of Magrini & Boschi, 2021)
prefac = -1i / (4*sqrt(2*pi)*vel^2);
arg = r * sqrt(omega.^2 / vel.^2 - 2*1i*alpha*omega/vel);
% Zero-order Hankel function of the 2nd kind
gf_f = prefac * besselh(0,2,arg);
gf_f(isnan(gf_f) | isinf(gf_f)) = 0; % deal with zero frequency...
gf_f = [gf_f, flip(conj(gf_f(2:end)))]; % build conjugate [0, +freq, -freq]

% Calculate Frequency-domain seismic response (eq 3 of Magrini & Boschi, 2021)
h = ones(size(gf_f)); % assume uniform strength as function of frequency
wavefield_f = h .* gf_f .* exp(1i*phi_rand);

% Get time-domain response
wavefield = real(ifft(wavefield_f));

% Plot example
if 0
    figure(999); clf;
    set(gcf,'position',[616   199   649   819]);
    
    t_arrival = r / vel;
    
    subplot(3,1,1); box on; hold on;
    plot([-1*flip(freq(2:end)), freq],abs(fftshift(gf_f)),'b','linewidth',2);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$G(x,x_S,\omega) = -\frac{i}{4\sqrt{2\pi}c^2} H_0^{(2)} \left(r \sqrt{\frac{\omega^2}{c^2} - \frac{2i\alpha\omega}{c}} \right)$','Interpreter','latex','fontsize',18)
    xlabel('Frequency (Hz)');
    
    subplot(3,1,2); box on; hold on;
    plot(t,real(ifft(gf_f)),'b','linewidth',2);
    plot(t_arrival*[1 1],[min(wavefield) max(wavefield)],'-r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$G(x,x_S,t)$','Interpreter','latex','fontsize',18)
    
    subplot(3,1,3); box on; hold on;
    plot(t,wavefield,'b','linewidth',2);
    plot(t_arrival*[1 1],[min(wavefield) max(wavefield)],'-r','linewidth',1);
    set(gca,'fontsize',15,'linewidth',1.5);
    title('$S(x,t) = F^{-1} \left[h(\omega) \, G(x,x_S,t) \, e^{i\phi} \right]$','Interpreter','latex','fontsize',18)
    xlabel('Time (sec)');
    
%     pause
    drawnow
end

end

