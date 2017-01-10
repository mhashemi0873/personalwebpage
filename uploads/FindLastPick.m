function LP = FindLastPick(t, f)
L=numel(f);
isFind=0;
for i=1:L
    if(t < f(i))
        LP = f(i-1);
        isFind=1;
        break;
    end;
end;
if(isFind == 0)
    LP = f(L);
end;
