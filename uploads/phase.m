%FitzHugh-Nagumo Euler version
clear all
clc
dt=.01;fn=121;in=0;
tspan=in:dt:fn;
endstep=((fn-in)/dt);
I=0.5;
a=.7;
b=.8;
eps=.08;
v=[];
u=[];
f=[];
ff=[];
v(1) = 0;
u(1) = 0;

for i=1:endstep
    v(i+1) = v(i) + dt*(v(i) - (v(i)^3/3) - u(i) + I)  ;
    u(i+1) = u(i) + dt*(eps*(v(i) + a - b*u(i)));
end

        f=Spike(v, endstep);
        p=numel(f);
                     
                   if(p <3)
                   fer=0;
                   else
                   m=(f(p-1)-f(1))/(p-2);
                   fer=1/(dt*m);
                   end
                   T=1/fer;
                 
                   counter=1;
                   for j=f(1):endstep
                       phi(counter)=(2*pi*dt)/T * (j - FindLastPick(j, f));
                       counter=counter+1;
                   end;
                   
                   x=[f(1):endstep];
%                    plot(v,u)
                    subplot(2,1,1)
                    plot(tspan,v)
                    subplot(2,1,2)
                    plot(dt*x,phi);
                   