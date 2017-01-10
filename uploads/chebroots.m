clc
clear
%dX/dt=Ax(t)+BX(t-tau)

B=-1;A=.5;tau=1;


N=100; n=length(A); % Discretization nodes N and size of DDE n
D=-cheb(N-1)*2/tau;
S=[kron(D(1:N-1,:),eye(n));[B,zeros(n,(N-2)*n), A]];
s=eig(S)

hold on
plot (real(s),imag(s),'b*')
grid on
xlabel('REAL')
ylabel('IMAGINARY')
axis([-10 0 -400 400 ])
title('dX/dt=AX(t)+BX(t-\tau)     B=-1,A=0.5,\tau=1')