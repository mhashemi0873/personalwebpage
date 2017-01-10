clc
clear

gNa = 120; gK = 36; gL = .3;
vNa = 115; vK = - 12; vL =  10.6; % neocortical pyramidal neurons 
current=7;

in=-11;fn=100;dt=.1;
endstep=((fn-in)/dt);

x=in:dt:fn;
y=[];
alpha_n=[];
alpha_m=[];

beta_n=[];
beta_m=[];
minf=[];
                                          
                                          alpha_m=@(x)  (.1*(25-x))/((exp((25-x)/10)-1));
                                          beta_m=@ (x)  4*(exp(-x/18));
                                       
                                          minf=@(x)      alpha_m(x)/(alpha_m(x)+beta_m(x));
                                          

f=@(y) -gNa*(minf(in)^3)*(0.89-1.1*y)*(in-vNa)-gK*(y^4)*(in-vK)-gL*(in-vL)+current;


y(1)=fsolve(f,0);

for i=2:endstep+1
    
                                         dx = in + dt*(i-1);
                                         f=@(y) -gNa*(minf(dx)^3)*(0.89-1.1*y)*(dx-vNa)-gK*(y^4)*(dx-vK)-gL*(dx-vL)+current;

                                         y(i)=fsolve(f,y(i-1));
end

for  i=1:endstep+1


                                        
                                          alpha_n(i)= (0.01*(10-x(i)))/((exp((10-x(i))/10)-1));
                                          beta_n(i)=  0.125*(exp(-x(i)/80));
                                          z(i)=alpha_n(i)/(alpha_n(i)+beta_n(i));
end


plot(x,y,'--')
hold on

plot(x,z,'r--')

%axis([-2.5 2.5 -2 2.5])  
xlabel('v');
ylabel('n');
hold on