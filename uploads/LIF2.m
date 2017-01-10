clc
clear

b=1.2;  % neocortical pyramidal neurons
dt=.01;
in=0;fn=20;
endstep=((fn-in)/dt);
tspan = in:dt:fn;
ref=50;
v=[];
v(1)=0; 
           
    for i=1:endstep
                  
                  v(i+1)=v(i)+(dt)*(-v(i)+b);
                  if v(i+1)>=1.1
                    v(i+1)=0;                                        
                  end;
     end;
   
     
 [pks,locs] = findpeaks(v);
     
 for i=1:numel(locs)-1
     
     
     for j=locs(i)+ref:locs(i+1)
                                     v(j+1)=v(j)+(dt)*(-v(j)+b);
                                      if v(j)>1 ; 
                                      v(j)=0 ;
                                      end;
     end
     
     for j=locs(i+1):locs(i+1)+ref
             v(j)=0;
     end
 end

t=[locs(1):locs(end)+locs(1)];

hold on
vv=v(locs(1):end);
tt=[0:length(vv)-1];
plot(dt*tt,vv)