function [X_turb,Y_turb,Z_turb,u] = turbulence_generator(filename,n1,n2,n3,Lx,Ly,Lz,H,umean)

% close all
% clear all
% This is an example of how to load in a Mann box file and plot the contour and corresponding PSD. 
% Use with care! I make many errors!

%Mann box file: % 'Turbulence_generator/sim1.bin'
fid=fopen(filename);

% %Number of points in box
% n1=4096; % z-direction
% n2=32; % y-direction
% n3=32; % x-direction
% 
% %Size of box (make sure to change to your setup!)
% Lz=3685.5;
% Ly=180;
% Lx=180;

deltay=Ly/(n2-1);
deltaz=Lz/(n1-1);
deltax=Lx/(n3-1);

% Time step
deltat = deltaz/umean;
% disp(deltat)

%Create matrix
uraw=fread(fid,'single');
itael=0;
u = zeros(n1,n2,n3);
for i=1:n1
    for j=1:n2
        for k=1:n3
            itael=itael+1;
            u(i,j,k)=uraw(itael);
        end
    end
end

%% Create plane coords (remember to shift plane to rotor position in your BEM code!)
X_turb = (0:n3-1)*deltax + (H-Lx/2);
Y_turb = (0:n2-1)*deltay - (Ly/2);
Z_turb = (0:n1-1)*deltaz;

end

% %Time vector
% time=deltat:deltat:n1*deltat;
% 
% %Signal to process. (Just picked a random point on plane)
% sig=u(:,20,20);
% 
% % PSD.
% t_start_PSD = 0;
% t_end_PSD = Inf;
% k_start_PSD = find(time>=t_start_PSD, 1, 'first')
% k_end_PSD   = find(time<=t_end_PSD,   1, 'last')
% timesim_PSD = time(k_start_PSD:k_end_PSD);
% sampling_frequency = 1/deltat;
% nyquist_frequency = sampling_frequency/2;
% n_parts_PSD = 2;
% window = hann(round(length(timesim_PSD)/n_parts_PSD), 'periodic');
% noverlap = round(length(window)/2);
% nperseg = length(window);
% nfft = nperseg;
% % Apply Welch's method.
% [sig_PSD, frequency] = pwelch(sig(k_start_PSD:k_end_PSD), window,  ...
%     noverlap, nfft, sampling_frequency, 'psd');
% % Plot.
% figure();
% loglog(frequency, sig_PSD, 'LineWidth', 2, 'DisplayName', 'Blade 1');
% % xlim([0 5]);
% xlabel('frequency [Hz]')
% ylabel('PSD(u)')
% 
% %     