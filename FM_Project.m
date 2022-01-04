%% cleaning
clc;
clear variables;
close all;

%% Reading Audio
% read the audio file 
[message , FS] = audioread ('eric.wav') ;

%sound time domain plot
timeAxis = linspace(0, length(message)/FS , length(message)) ;
figure
subplot(2,1,1)
plot (timeAxis , message,'r') ; xlabel('time(s)') ; title('Message Time Domain Waveform');
%message spectrum plot
Fdom1 = linspace(-FS/2 , FS/2 , length(message));
message_freq_domain = fftshift(fft(message)) ;
message_mag = abs(message_freq_domain) ;
subplot(2,1,2); 
plot (Fdom1, message_mag,'r') ; xlabel('frequency(Hz)') ; title('Message Spectrum');

%% Filtering
n = FS/(length(message)-1) ;
% first interval
freq1 = -FS/2 : n : -4000 ;
mag1 = zeros(1, length(freq1)) ;
% second interval
freq2 = -4000: n : 4000 ;
mag2 = ones(1, length(freq2)) ;
% third interval
freq3 = 4000 : n : FS/2 ;
mag3 = zeros(1, length(freq3));
% filtering 
filter = [mag1 mag2 mag3];
filtered_message_F = filter' .* message_mag ;

% converting from frequency domain to time domain 
size(filtered_message_F); size(message) ;
filtered_message_T = real(ifft(ifftshift(message_freq_domain))) ;

% PLOTTING
figure
subplot(2,1,2);
plot(Fdom1 , filtered_message_F,'g') ; xlabel('frequency(Hz)') ; title('filtered message spectrum'); 
subplot(2,1,1);
plot(timeAxis , filtered_message_T,'g') ; xlabel('time(s)')  ; title('filtered message time domain') ;

%play sound
%sound(filtered_message_T , FS) ;

%% NBFM Modulation
% Modulation equation
%   s(t) = A cos(2 pi Fc t) - A Kf sin(2 pi Fc t) integration(message)

Fc = 100000;                 % Carrier Frequency
FS2 = 5 * Fc;                 % Sampling Frequency
Kf = 2 * pi * 0.01;           % Frequency Sensitivity
A = 10;                           % Carrier Amplitude

% resampling the message
message_resampled = resample(filtered_message_T , FS2 , FS);

% generating linspaces(defining time and freq intervals)
t = linspace (0 ,  (length(message_resampled)) / FS2 , length(message_resampled)) ;
Fdom2 = linspace( (-FS2)/2 , FS2/2 , length(message_resampled) ) ;

% generating the signal
NBFM_signal = A * cos(2*pi*Fc* t)' - A * Kf * sin(2*pi*Fc* t)' .* cumsum(message_resampled);
%converting from time to freq. domain
NBFM_signal_spectrum = fftshift ( fft (NBFM_signal)) ;

% PLOTTING
figure
subplot(2,1,1);
plot ( t,NBFM_signal,'m') ; xlabel('time(s)') ;  title('NBFM message in time domain');
subplot(2,1,2);
plot(Fdom2, abs(NBFM_signal_spectrum),'m') ; xlabel('frequency(Hz)') ; title('NBFM magnitude');

%% Demodulation
% demodulation is achieved by using a differentiator followed by 
% an envelope detector

% envelope detector
out2 = abs(hilbert(NBFM_signal));
out2_f = fftshift(fft(out2)) ;
%plot
figure
subplot(2,1,1)
plot(t, out2, 'c') ; xlabel('time(s)') ; title(' filtered envelope signal in time domain ') ;
subplot(2,1,2)
plot(Fdom2, abs(out2_f), 'c') ; xlabel('freq(Hz)') ; title('filtered envelope signal in freq domain') ;

% differentiator
out1 = [0; diff(out2)];
%plot
figure
plot(t, out1); xlabel('time(s)') ;title('Demodulated Message') ;
ylim([-0.5 0.5]);
%playing final sound
sound(3 * resample(out1 , FS , FS2), FS)