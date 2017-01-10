clc
clear all
format long

D = 0.1;
gama = 1;
in=0;fn=10;
dt = 0.001;
endstep =(fn-in)/dt;
tspan=in:dt:fn;

v0 = 0;
realization = 100;
L= 10000;
Fs = 20;
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2);

v = zeros(realization,L);
v(:,1) = v0;


for i=1:realization
    
    noise = wgn(1,endstep+1,2*D);
    
    for j=1:endstep
        v(i,j+1) = v(i,j)*(1-gama*dt)+sqrt(2*D*dt)*noise(j); 
    end  
   
    Fourier(i,:) = fft(v(i,:),NFFT)*dt;
    PS(i,:) = Fourier(i,:).*conj(Fourier(i,:));
    C(i,:) = ifft(PS(i,:));
    
end

PowerS = sum(PS)/(L*dt*realization);
Cor = sum(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1);
plot(tspan,v(1,:));
xlabel('Time');ylabel('v ( t )');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2);
plot(f(100:length(f)),300*Fourier(1,100:NFFT/2));
xlabel('Frequency');ylabel(' Fourier transform ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
plot(f(10:length(f)),Cor(10:NFFT/2),'o');
xlabel('Time delay');ylabel('Correlation function');
axis([0 10 -0.005 0.1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4);
plot(f(100:NFFT/2),PowerS(100:NFFT/2),'o');
xlabel('Frequency');ylabel('Power spectrum');
axis([0 10 0 0.0001]);
