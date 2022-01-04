clc
clear all
close all


%%
%%%PART 1_Reading Audio%%%

%read the audio file 
[signal , FS] = audioread ('eric.wav') ;

%sound time domain plot
timeAxis = linspace(0, length(signal)/FS , length(signal)) ;
figure
subplot(2,1,1)
plot (timeAxis , signal,'r'); xlabel('time(s)') ; title('Time Domain Waveform signal');
%signal spectrum plot
Fdom1 = linspace(-FS/2 , FS/2 , length(signal));
signal_freq_domain = fftshift(fft(signal)) ;
signal_mag = abs(signal_freq_domain) ;
subplot(2,1,2); 
plot (Fdom1, signal_mag,'r') ; title('frequency signal');

%play the audio file
%sound (signal , Fs) ;


%%
%_____________________________
%FILTER 

n = FS/(length(signal)-1) ;
%first interval
freq1 = -FS/2 : n : -4000 ;
mag1 = zeros(1, length(freq1)) ;
%second interval
freq2 = -4000: n : 4000 ;
mag2 = ones(1, length(freq2)) ;
%third interval
freq3 = 4000 : n : FS/2 ;
mag3 = zeros(1, length(freq3));
%filtering 
filter = [mag1 mag2 mag3];
filtered_signal_F = filter' .* signal_mag ;
%converting from frequency domain to time domain 
size(filtered_signal_F); size(signal) ;
filtered_signal_T = real(ifft(ifftshift(signal_freq_domain))) ;

%PLOTTING
figure
subplot(2,1,2);
plot(Fdom1 , filtered_signal_F,'g') ; title('filtered signal spectrum'); 
subplot(2,1,1);
plot(timeAxis , filtered_signal_T,'g'); title('filtered signal time domain') ;
%playing audio
%sound(filtered_signal_T , FS)

%%
%________________________________
%Modulation
fc = 100000 ;
FS2 =5*fc ;
kf = 0.02 ;

t = linspace (0 ,  (length(filtered_signal_T)) / FS , length(filtered_signal_T)) ;
FM = 10*cos(2*pi*fc*t' +2*pi*kf*cumsum(filtered_signal_T)) ;

Fdom2 = linspace(-FS2/2 , FS2/2 , length(FM)) ;
%converting from time into freq domain
FM_spec = fftshift(fft(FM));
FM_mag = abs(FM_spec);

%plot
figure
subplot(2,1,1) ;
plot(Fdom2 , FM_mag); title('NBFM spectrum');
subplot(2,1,2) ;
plot(t , FM); title('NBFM time domain');

%%
%________________________________
%Demodulation
envelope = FM ;
X = resample(envelope, length(filtered_signal_T)/FS , length(envelope)) ; 

Y= diff(envelope);
envelope = abs(hilbert(Y)) ;

figure 
plot (envelope) ;

sound(envelope , FS)

%  envelope = abs(hilbert(FM)) ;
%  envelopeF = fftshift(fft(envelope)) ;
 
 % sig_lpf = (4000*length(envelopeF))/ FS ;
 % envelopeF([ 1 : floor(length(envelopeF)/2-sig_lpf) ceil(length(envelopeF)/2 +sig_lpf):end]) =0 ;
 
 

% e = abs(hilbert(s1));
% E=fftshift(fft(e)); %turn it into frequency domain
% lpfe = 4000 * length(E) / fs;
% E([ 1:floor(length(E)/2 - lpfe) ceil(length(E)/2 + lpfe):end ]) = 0;
% re = real(ifft(ifftshift( E)));
% sig = diff(re);
% sound(sig,fs) % to play audio with fs frquenc



%mod_sig2 = resample(smt, length(signal) , length(smt)) ;

% diff_smt = diff(smt) ;
% env_det = real ( hilbert (diff_smt) ) ;
% 
% figure
% plot (env_det) ; title('demodulated signal') ;
%sound(env_det);



%%
%%% PART 2_Filter %%%
% %defining the filter interval
% a = ((length(signal_freq_mag))/Fs) * (Fs/2 -4000) ;
% b = ((length(signal_freq_mag))/Fs) * (Fs/2 +4000) ;
% % rounding elements to the nearest integer for more accuracy results
% c = round(a) ;
% d = round(b) ;
% signal_freq_mag([1:c d:length(signal_freq_mag)]) = 0 ;
% new_signal_mag = abs(signal_freq_mag) ;
% %plot filtered signal spectrum
% figure
% subplot(2,1,1);
% plot (Fdomain , new_signal_mag,'g'); title('Filtered Signal Spectrum'),ylabel('amplitude');
% 
% %transferring from frequency to time domain
% new_signal_time_domain = real(ifft(ifftshift(signal_freq_domain) ) ) ;
% 
% %plot filtered signal in time domain
% subplot(2,1,2);
% plot(timeAxis,new_signal_time_domain,'g') ; title ('Filtered Signal Time Domain'); xlabel('time(s)');
% 
% %play the audio file
% %sound(new_signal_time_domain,Fs);
%%
%ev_det = s1 .sin(wc(t2)');
% fc = 100000 ;               %Carrier Frequency
% FS2 = 5*fc ;                     %Sample Frequency
% kf = 0.02;                        %Frequency Sensitivity
% %define time interval & frequency interval
% new_S = resample(filtered_signal_T , FS2 , FS);
% t = linspace (0 ,  (length(new_S)) / FS2 , length(new_S)) ;
% Fdom2 = linspace( (-FS2)/2 , FS2/2 , length(new_S) ) ;
% %FM modulation
% FM = 10*cos(2*pi*fc*t' +(2*pi*kf*cumsum(new_S) ) ) ;
% FM_spectrum = fftshift ( fft (FM)) ;
% %figure 
% %plot(t,trans_signal);
% % Fm_spectrum = abs ( fftshift ( fft (Fm) ) ) ;
% figure
% subplot(2,1,1);
% plot ( t,FM); title('FM modulation signal in time domain');
% subplot(2,1,2);
% plot(Fdom2, abs(FM_spectrum)); title('NBFM magnitude');
