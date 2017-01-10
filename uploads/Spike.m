function f=Spike(v, endstep)
f=[];
p=1;
k=1;
flag = 0;
if v(k+1) > v(k)
    flag = 1;
end;

for k=2:endstep-1
      if flag == 1
          if v(k+1) > v(k)
              continue;
          else
              f(p)=k;
              p=p+1;
              flag = 0;
          end;
      else
          if v(k+1) < v(k)
              continue;
          else
              flag = 1;
          end;
      end;                           
end;