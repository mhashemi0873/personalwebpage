clc
clear
%dx/dt=ax(t)+bx(t-tau)
 
b=-1;a=.5;tau=1;

                    
  for k=-50:50;
  s=((1/tau)*lambertw(k,b*tau*exp(-a*tau)))+a;
  hold on
  plot(s,'go','markersize',8)
  end
                   
grid on
xlabel('REAL')
ylabel('IMAGINARY')
axis([-10 0 -400 400 ])
title('x^.=ax+bx(t-\tau)     b=-1,a=0.5,\tau=1')