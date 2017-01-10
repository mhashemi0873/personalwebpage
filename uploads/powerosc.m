clc
clear 

D = .5;
gama =.5;
omega=2*pi*1;
realization = 10;
Fs = 1000;                         % Sampling frequency
dt = 1/Fs;                         % Sample time
L = 20000;                        % Length of signal
              
NFFT = 2^nextpow2(L);             % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);
v10 =10;                         %initial condition
v20 =10;                         %initial condition

v1= zeros(realization,L);
v2= zeros(realization,L);

v1(:,1) = v10;
v2(:,1) = v20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:realization
    
    
    for j=1:L
        
        v1(i,j+1)=  v1(i,j) + dt*v2(i,j); 
        v2(i,j+1) = v2(i,j) + dt*(-omega^2*v1(i,j)-gama*v2(i,j))+sqrt(2*D*dt)*randn; 
    end
    
    
    M=length(v1);
    
    H = 1.05*(1 - cos(2*pi*(1:M)'/(M+1)));    %
    v0(i,:) = v1(i,:).* H' ; 
  
    Y(i,:) = abs(fft(v0(i,:),NFFT)).^2 *(dt/NFFT);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=0;fn=20;

PowerS = mean(Y);
z=PowerS(1:NFFT/2+1);

Hs = spectrum.periodogram('Hamming');
psd(Hs,mean(v1),'Fs',Fs,'NFFT',1024,'SpectrumType','onesided')
Hpsd = dspdata.psd(z,'Fs',Fs);

subplot(2,2,1)
hold on
plot(Hpsd);
hold on
plot(f,10*log10(z),'r--')
%set(gca, 'yscale', 'log')
xlim([0 fn])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=0:dt:fn;
w=w*2*pi;
freq=0:dt:fn;

Den=(omega^2-w.^2)+1i*gama*w;
Re=real(Den);Im=imag(Den);
Ginv=1./Den;

P=(2*D)*(1./((((w.^2)-(omega^2)).^2)+((gama.*w).^2)));
%P=(2*D)*Ginv.*conj(Ginv);

plot(freq,10*log10(P),'g');
%set(gca, 'yscale', 'log')
ylabel('Power spectrum (a.u)','interpreter','latex','fontsize',14);
xlabel('Frequency (Hz)','interpreter','latex','fontsize',14);
set(gca, 'fontsize',14);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J=[0  1;-omega^2 -gama];
lambda=real(eig(J))+1i*imag(eig(J))/(2*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ww=[(1i*gama/2)+sqrt(omega.^2-((gama.^2)./4))/(2*pi) 
       (1i*gama/2)-sqrt(omega.^2-((gama.^2)./4))/(2*pi)]
   
lambda=1i.*ww

f_max=sqrt(omega.^2-((gama.^2)./2))/(2*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
time=0:dt:L*dt;
hold on
plot(time,mean(v1));
ylabel('Amplitude (a.u)','interpreter','latex','fontsize',14);
xlabel('Time [sec]','interpreter','latex','fontsize',14);
set(gca, 'fontsize',12);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
 [S,F,T]=spectrogram(mean(v1), Hann(15000),7000,2000,Fs,'yaxis');  
%[S,F,T]= spectrogram(mean(v1),blackman(100000), 25000, 20000,Fs,'yaxis');                                    
                                            %Compute fft over 1000 pts.
S = abs(S);
bmin=max(max(abs(S)))/100000000;
hold on
imagesc(T,F,10*log10(max(abs(S),bmin)/bmin));
%imagesc(T,F,10*log10(S/min(S(:)))); 

%  surf(T,F,10*log10(max(abs(S),bmin)/bmin),'edgecolor','none')
%  axis tight;
%  view(0,90);

colorbar;
colormap(jet(256))
axis xy
ylim([0 10])
ylabel('Frequency [Hz]','interpreter','latex','fontsize',14);
xlabel('Time [sec]','interpreter','latex','fontsize',14);
set(gca, 'fontsize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt=0:dt:fn;

Cxx = xcorr(mean(v1),fn/dt);
[junk,index] = max(Cxx);

subplot(2,2,4);
hold on
plot(tt,Cxx(index:end)/Cxx(index),'r')

omega_0=sqrt(omega.^2-((gama.^2)./4));
CC=(D./(gama*(omega^2))).*exp(-(gama/2).*tt).*(cos(omega_0.*tt)+((gama/(2*omega_0))*sin(omega_0.*tt)));

hold on
plot(tt,CC/CC(1),'b')
xlim([0 10])
set(gca, 'fontsize',14);
hold on
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [ACF, Lags, Bounds] = autocorr(mean(v1),fn/dt);
% subplot(2,1,2);
% plot(tt,ACF,'go');

% Rx=autocorrelation(mean(v1));
% plot(tt,Rx/Rx(1),'ro');
% grid;
% subplot(2,1,2)
