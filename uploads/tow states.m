clc
clear all
format long

T = [10,20,30];     % Simulation Time
dt = 0.01;          % Time Step
tau=[0.1,.5,1];   % Refractory Period

rup = 0.5/dt;       % For leaving the Up State
rdn = 1-rup;        % For leaving the Down State

sigmaup = 1;         % Up State
sigmadn = -1;        % Down State



 

sigma(1) = sigmaup;  % Initial Condition
time(1) = 0;
sig = sigmaup;
       
       
np=1;

  for counter=1:size(tau,2)  
      
      jj=1;ii=1;            % Loop Counter
      
nt = T(counter)/tau(counter);        % Number of Steps
 r = rand(1,nt+1);    % Uniform Distribution Random Number
 
 TD=tau(counter);
 
  for i=0:TD:T(counter)-TD
    
           for j=i:dt:i+TD
               
           sigma(ii) = sig;
           time(ii)=j;
           ii=ii+1;
           end
                  if r(1,jj)<rup*dt
                      sigma(ii) = sigmaup;
                      sig = sigmaup;
                  else
                      sigma(ii) = sigmadn;
                      sig = sigmadn;
                  end
      jj=jj+1;
  end
     
     time(ii)=T(counter);


     subplot(3,1,np);
     plot(time,sigma);axis([0 T(counter) sigmadn-1 sigmaup+1])
     ylabel(['\tau_{D}=', num2str(TD)]);
     if ii==3
         xlabel('Time');
     end
     np=np+1;
   end