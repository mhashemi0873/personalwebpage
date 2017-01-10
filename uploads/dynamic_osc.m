clc
clear 

D = .5;
gama =.5;
omega=2*pi*1;
realization = 10;
Fs = 1000;                     % Sampling frequency
dt = 1/Fs;                     % Sample time
L = 100000;                     % Length of signal
              
NFFT = 2^nextpow2(L);           % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);
v10 =.1;                         %initial condition
v20 =.1;                         %initial condition

v1= zeros(realization,L);
v2= zeros(realization,L);

v1(:,1) = v10;
v2(:,1) = v10;


for i=1:realization
    
    
    for j=1:L
        
        v1(i,j+1)=  v1(i,j) + dt*v2(i,j); 
        v2(i,j+1) = v2(i,j) + dt*(-omega^2*v1(i,j)-gama*v2(i,j))+sqrt(2*D*dt)*randn; 
    end
    
    
    M=length(v1);
    H = (1 - cos(2*pi*(1:M)'/(M+1)));    
    %H = hanning(M);
 
    v0(i,:) = v1(i,:).*H' ; 
  
    Y(i,:) = abs(fft(v0(i,:),NFFT)).^2 *(dt/NFFT);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=0;fn=10;

PowerS = mean(Y);
z=PowerS(1:NFFT/2+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
hold on
plot(f,z,'r')
%set(gca, 'yscale', 'log')
xlim([0 fn])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=0:dt:fn;
w=w*2*pi;
freq=0:dt:fn;

Den=(omega^2-w.^2)+1i*gama*w;
Re=real(Den);Im=imag(Den);
Ginv=1./Den;

P=(2*D)*(1./((((w.^2)-(omega^2)).^2)+((gama.*w).^2)));
%P=(2*D)*Ginv.*conj(Ginv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
plot(freq,P);
set(gca, 'yscale', 'log')
xlabel('Frequency (Hz)','interpreter','latex','fontsize',12);
ylabel('Power spectrum (a.u.)','interpreter','latex','fontsize',12);
set(gca, 'fontsize',14);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[0  1;-omega^2 -gama];

[V,D]=eig(A);

lambda=real(eig(A))+1i*imag(eig(A));

lambda1=real(D(1,1))+1i*imag(D(1,1))/(2*pi)
lambda2=real(D(2,2))+1i*imag(D(2,2))/(2*pi)

V1=V(:,1)
V2=V(:,2)

J=[0-lambda(1)  1;-omega^2 -gama-lambda(1)];
DET=det(J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ww=[(1i*gama/2)+sqrt(omega.^2-((gama.^2)./4))/(2*pi) 
       (1i*gama/2)-sqrt(omega.^2-((gama.^2)./4))/(2*pi)]
   
lambda=1i.*ww

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_max=sqrt(omega.^2-((gama.^2)./2))/(2*pi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C=V\[v10 v20]';

t=0:dt:60;

X1=C(1)*V(1,1).*exp(lambda1.*t)+C(2)*conj(V(1,1).*exp(lambda1.*t));
X2=C(1)*V(2,1).*exp(lambda1.*t)+C(2)*conj(V(2,1).*exp(lambda1.*t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
hold on
plot(X1,X2)
xlabel('X(t)','interpreter','latex','fontsize',14);
ylabel('$$X^{\cdot}(t)$$','interpreter','latex','fontsize',14);
set(gca, 'fontsize',14);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
hold on
plot(t,X1)
hold on
plot(t,X2)
xlabel('Time','interpreter','latex','fontsize',14);
ylabel('$$ X(t)~ \&  ~X^{\cdot}(t) $$','interpreter','latex','fontsize',14);
set(gca, 'fontsize',14);
box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;
b=(-omega^2)./(lambda1-(gama));

%a=V1(1,1);
%b=V1(2,1);

cof=1./sqrt((a.^2)+(b.^2));

u=[a b]/norm([a b]);

R_l=real(lambda1);I_l=imag(lambda1);
%R_l=real(lambda2);I_l=imag(lambda2);

R_1=real(u(1)); I_1=imag(u(1));
R_2=real(u(2)); I_2=imag(u(2));


syms t

F1=@ (t) exp(2*R_l.*t).*((cos(I_l.*t)).^2);
F2=@ (t) exp(2*R_l.*t).*((sin(I_l.*t)).^2);
F3=@ (t) exp(2*R_l.*t).*(cos(I_l.*t).*sin(I_l*t));

if I_l==0
I1=int(F1,t,0,(2*pi));
I2=int(F2,t,0,(2*pi));
I3=int(F3,t,0,(2*pi));
else
I1=int(F1,t,0,(2*pi)/I_l);
I2=int(F2,t,0,(2*pi)/I_l);
I3=int(F3,t,0,(2*pi)/I_l);
end

I1=eval(I1);
I2=eval(I2);
I3=eval(I3);

x_1= 4*(((R_1^2)*I1)+((I_1^2)*I2)-(2*(R_1*I_1)*I3));
x_2= 4*(((R_2^2)*I1)+((I_2^2)*I2)-(2*(R_2*I_2)*I3));



X_n=[x_1 x_2]'

W1=x_1/sum(X_n);
W2=x_2/sum(X_n);


Y1=W1
Y2=W2

hold off
subplot(2,2,4)

Population = {'V_1', 'V_2','interpreter','latex'};
Y = 100*[Y1 
         Y2];
 
bar(Y)
colormap (jet(11))
grid on
set(gca, 'XTickLabel', Population,'fontsize',14);
xlabel('','interpreter','latex');

ylabel('$$ Contribution ~in ~power \% $$','interpreter','latex','fontsize',12);

print -dpdf plot.pdf
%print -depsc plot.eps