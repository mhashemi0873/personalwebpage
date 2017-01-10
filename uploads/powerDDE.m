clc
clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  dx/dt=ax(t)+bx(t=tau)  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%stability conditions:
%%% if b>0 then  a<-b
%%% if b<0 then a<-b*cos(omega*tau) : omega=-bsin(omega*tau)%%  0<omega< pi/tau
%%%   or
%%%  a<-b and tau*sqrt(b^2-a^2)< acos(-a/b)

%%%% omega^2=b^2-a^2 %%%%%%%

D = .1;
omega=5*(2*pi);
tau=.05;
epsilon=1;
a=(omega/tan(omega*tau))-epsilon;
b=-omega./sin(omega*tau);

realization = 50;
Fs = 1000;                       % Sampling frequency
dt = 1/Fs;                       % Sample time
L = 100000;                      % Length of signal
NFFT = 2^nextpow2(L);            % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);
v0 = 1;
v= zeros(realization,L);
v(:,1) = v0;
delay=tau/dt;

for i=1:realization
    
        noise = wgn(1,L+1,1);
        
    for j=1:delay+1
        v(i,j)=v0;
    end
    for j=delay+1:L
         
        v(i,j+1) = v(i,j) + dt*(a*v(i,j)+b*v(i,j-delay))+sqrt(2*D*dt)*randn; 
    end 
    
    M=length(v);
    H1 = (1 - cos(2*pi*(1:M)'/(M+1)));    %
    H2 = taylorwin(M);
    H=(H1+H2)/2;
    v(i,:) = v(i,:).* H' ; 
  
    Y(i,:) = abs(fft(v(i,:),NFFT)).^2 *(dt/NFFT);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211);
plot(v(1,:));
xlabel('Time (msec)');ylabel('v ( t )');
ylabel(' v ( t )','interpreter','latex','fontsize',18)
xlabel(' Time  ','interpreter','latex','fontsize',18)
set(gca, 'fontsize',16);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(212);
PowerS = sum(Y)/realization;
in=0;fn=50;
z=10*log10(PowerS(1:NFFT/2+1));
plot(f,z,'r')
xlim([0 fn])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
w=0:dt:fn;
w=w*2*pi;
P=(2*D)*(1./(((a+b*cos(tau*w)).^2)+((w+b*sin(tau*w)).^2)));
%or P=((2*(D.^2)/sqrt(2*pi)))*(1./((a^2)+(b^2)+(w.^2)+(2*a*b*cos(tau.*w))+(2*b*w.*sin(tau.*w))));
freq=0:dt:fn;
plot(freq,10*log10(P));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ylabel(' Spectral power (a.u.)','interpreter','latex','fontsize',14)
xlabel(' Frequency (Hz)  ','interpreter','latex','fontsize',16)
set(gca, 'fontsize',16);
box on

