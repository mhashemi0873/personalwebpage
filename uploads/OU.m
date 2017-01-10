clc
clear all
format long

D=0.1;
gama =1;
in=0;fn=5;
dt = 0.01;
endstep =(fn-in)/dt;
tspan=in:dt:fn;
realization = 10000;
v0 = 1;%initial velocity
v = zeros(realization,endstep);
v(:,1) = v0;

for i=1:realization
    
    noise = wgn(1,endstep+1,2*D);
    
    for j=1:endstep
        v(i,j+1) = v(i,j)*(1-gama*dt)+sqrt(2*D*dt)*noise(j); 
    end
    
    if (rem(i,10000)==1)
       plot(tspan,v(i,:))
       hold on
    end
end