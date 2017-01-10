  
clear all
clc
dt=.1;fn=100;in=0;
tspan=in:dt:fn;
endstep=((fn-in)/dt);

I=.459;
a=.7;
b=.8;
eps=.08;
vv=[];
uu=[];
vv(1)=0;
uu(1)=0;
v=[];
u=[];
f=[];
ff=[];
delta=[];
phi0=[];
v(1) = 0;
u(1) = 0;
   


for i=1:endstep 
      v(i+1) = v(i) + dt*(v(i) - (v(i)^3/3) - u(i) + I )  ;
      u(i+1) = u(i) + dt*(eps*(v(i) + a - b*u(i)));
end

        f=Spike(v, endstep);
        p=numel(f);
                     
        if(p < 3)
        fer=0;
        else
        m=(f(p-1)-f(1))/(p-2);
        fer=1/(dt*m);
        end
        T=1/fer;
        temp = T/dt;
        disp(temp);
        
jj=1;
for T_0=1:1:temp-.1
    
    First = f(1)+T_0;
    L= First+5;   
    for ii=1:endstep 
                if(ii>First && ii< L)
                    Current=1.;                 
                else
                    Current = 0;
                end;
          vv(ii+1) = vv(ii) + dt*(vv(ii) - (vv(ii)^3/3) - uu(ii) + I+Current )  ;
          uu(ii+1) = uu(ii) + dt*(eps*(vv(ii) + a - b*uu(ii)));  
    end

    ff=Spike(vv, endstep);

    phi0(jj)= (2*pi*dt)/T * (First- FindLastPick(First, f));
    delta(jj)=0.01*(((f(end)-ff(end))/T));
    jj=jj+1;  
end;                
                     

%plot(v)
hold on
%plot(vv,'go')
%plot(phi);
%plot(phi,'r');
subplot(2,1,2)
 plot(phi0,delta,'b');
    
 