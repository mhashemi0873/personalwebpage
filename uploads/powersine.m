clc
clear
Fs = 1000;                   % Sampling frequency.
dt = 1/Fs;                   % Sample time (here is 1 msec).
L = 10000;                   % Length of signal (here is 10000 msec =10 sec).
t = (0:L)*dt;                % Time vector from zero to 10 sec.
NFFT = 2^nextpow2(L);        % Next power of 2 from length of y  (it is always fixed).
f = Fs/2*linspace(0,1,NFFT/2+1);   % the range of frquncy.
noise = 1.*randn(size(t));  %this niise. in below we add noise to systems.
v=10*sin(2*pi*5*t)+10*sin(2*pi*10*t)+noise;  % v is signal which should show two peaks at 5 Hz and 10Hz.

M=length(v);
    
H = (1 - cos(2*pi*(1:M)'/(M+1)));   % H = hanning(M);  %here we use a window.
vv = v.* H' ;
         
Y = abs(fft(vv,NFFT)).^2 *(dt/NFFT);  % power=abs(fft.^2) and (dt/NFFT) is normalization.
    
z=Y(1:NFFT/2+1);
PowerS=10*log10(z);       % here power is represented in 10 times of log in base of 10.



 

Hs = spectrum.periodogram('Hamming');            %this is the matlab function
psd(Hs,v,'Fs',Fs,'NFFT',1024,'SpectrumType','onesided')   
Hpsd = dspdata.psd(z,'Fs',Fs);

%%this is good for low frequency peaks Fs=100
%h = spectrum.periodogram('Hamming');
%Hpsd=psd(h,v,'Fs',Fs);

%%this is good for high frequency peaks Fs=10000
%h = spectrum.welch; 
%Hpsd=psd(h,y,'Fs',Fs); 



subplot(3,1,1)
hold on   
plot(v(1:1000))
set(gca, 'fontsize',14);
xlabel('t [ms]','interpreter','latex','fontsize',16);
ylabel('v [mv]','interpreter','latex','fontsize',16);
box on


subplot(3,1,2)
plot(Hpsd);
hold on
plot(f,PowerS,'r--')
%set(gca, 'yscale', 'log')
set(gca, 'fontsize',14);
xlabel('Frequency [Hz]','interpreter','latex','fontsize',14);
ylabel('Power [dB]','interpreter','latex','fontsize',16);
xlim([0 30])


subplot(3,1,3)
[S,F,T]=spectrogram(v,hamming(1000),100,400,Fs,'yaxis');  

                                           
S = abs(S);
bmin=max(max(abs(S)))/100000;
imagesc(T,F,10*log10(max(abs(S),bmin)/bmin));
%imagesc(T,F,10*log10(S/min(S(:)))); 

%  surf(T,F,10*log10(max(abs(S),bmin)/bmin),'edgecolor','none')
%  axis tight;
%  view(0,90);

colorbar;
colormap(jet(256))
axis xy
ylim([0 30])
ylabel('Frequency [Hz]','interpreter','latex','fontsize',12);
xlabel('Time [sec]','interpreter','latex','fontsize',16);
set(gca, 'fontsize',14);

