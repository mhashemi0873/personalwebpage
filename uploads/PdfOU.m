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
    
    if (rem(i,1000)==1)
       plot(tspan,v(i,:))
       hold on
    end
end

figure
n0=1;
t=[];
dv=.1;

for ii=1:endstep
    
    m(ii) = min(v(:,ii));
    M(ii) =max(v(:,ii));
    cell(ii) = abs(round((M(ii)-m(ii))/dv))+1;
    pdf = zeros(1,cell(ii)+1);
    
     for jj=1:realization
       
       wcell = round((v(jj,ii)/dv))-round(m(ii)/dv)+1;
       pdf(wcell)= pdf(wcell)+1;
      
     end
    
    x=m(ii):dv:M(ii);
    x1 = size(x);
    
    t(1) = 0.1/dt;t(2)=0.3/dt;t(3)=1/dt;t(4)=5/dt;
   
    
        if ii==t(n0)
        subplot(4,1,n0);
        plot(x,pdf(1:x1(2))/(realization*dv))
        axis([-1.5 1.5 0 2])

        c = dt*ii;
        title(['t = ',num2str(c)]);
        n0 = n0+1;
        end
    
end