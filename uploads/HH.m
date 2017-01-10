clc
clear

C=1;

gNa = 120; gK = 36; gL = .3;
vNa = 115; vK = - 12; vL =  10.6; % neocortical pyramidal neurons 

current=7;

dt=.01;

in=0;fn=100;

endstep=((fn-in)/dt);
tspan = in:dt:fn;



              v(1) = 0;
             Inoise(1)=0;
             
              alpha_n(1)=(0.01.*(10-0))/((exp((10-0)/10)-1));
              beta_n(1)=0.125.*(exp(-0/80));
                                          
              alpha_m(1)=(.1.*(25-0))/((exp((25-0)/10)-1));
              beta_m(1)=4.*(exp(-0/18));
                alpha_h(1)=0.07*(exp(0/20));
               beta_h(1)=1/((exp((30-0)/10)+1));                        
            
             
              n(1)=(alpha_n(1)/(alpha_n(1)+beta_n(1)));
              m(1)=(alpha_m(1)/(alpha_m(1)+beta_m(1)));
              h(1)=(alpha_h(1)/(alpha_h(1)+beta_h(1))); 
% 
%               alpha_n(1)=0;
%               alpha_m(1)=0;
%               beta_m(1)=0;
%               alpha_h(1)=0;
%               beta_h(1)=0;
% 
%               n(1)=0.;
%                m(1)=0.0;
%                h(1)=0.;
%     
D=20;
   noise=wgn(1,endstep+1,1); 
   
    for i=1:endstep        
        
        
                                          alpha_n(i)=(0.01*(10-v(i)))/((exp((10-v(i))/10)-1));
                                          beta_n(i)=0.125*(exp(-v(i)/80));
                                          
                                          alpha_m(i)=(.1*(25-v(i)))/((exp((25-v(i))/10)-1));
                                          beta_m(i)=4*(exp(-v(i)/18));
                                          
                                          alpha_h(i)=0.07*(exp(-v(i)/20));
                                          beta_h(i)=1/((exp((30-v(i))/10)+1));
                                          
                                          n(i+1)=n(i)+dt*(alpha_n(i)*(1-n(i))-beta_n(i)*n(i));
                                          m(i+1)=m(i)+dt*(alpha_m(i)*(1-m(i))-beta_m(i)*m(i));
                                          h(i+1)=h(i)+dt*(alpha_h(i)*(1-h(i))-beta_h(i)*h(i));
                                          
                                          Inoise(i+1)=Inoise(i)+(dt/2)*(-Inoise(i))+(sqrt(dt*(2*D)))*noise(i);
                                          
                                          INa(i)=(gNa)*(m(i)^3)*h(i)*(v(i)-vNa);
                                           IK(i)=(gK)*(n(i)^4)*(v(i)-vK);
                                           IL(i)=(gL)*(v(i)-vL);
                                      
                                      
                                           v(i+1)=v(i)+(dt/(C))*(-INa(i)-IK(i)-IL(i)+current+0*Inoise(i));
                                           
       
     end;


plot(tspan,v)

% f=Spike(v,endstep);
% p=numel(f);
% 
            p=1;
            for k=1:endstep-1
                   %if v(k)<-5 && v(k+1)>-5 && abs(v(k)-v(k+1))<.1
                   if v(k)*v(k+1)<0 && v(k)>0
                       f(p)=k;
                       p=p+1;
                   end 
            end
            if(p < 5)
            fer=0;
            else
            mm=(f(p-1)-f(1))/(p-2);
            fer=1/(dt.*mm);
            end
            T=1/fer