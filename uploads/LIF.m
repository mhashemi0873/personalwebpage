clc
clear

b=1.5;  % neocortical pyramidal neurons
dt=.01;
in=0;fn=20;
endstep=((fn-in)/dt);
tspan = in:dt:fn;
ref=50;
v=[];
v(1)=0; 
 ref=50;          
flag = 0;

    for i=1:endstep
        
              if(flag > i)
                  v(i+1) = 0; 
              else
                  
                  v(i+1)=v(i)+(dt)*(-v(i)+b);
                  if v(i+1)>=1.01
                    v(i+1)=0; 
                    flag = i+ref;                                        
                  end;
              end;
    
     end;
 plot(tspan,v)    
     

