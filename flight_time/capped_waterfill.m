function p = capped_waterfill(a,cap,Pmax,nu)
% Fixed-assignment KKT solution: one interference cap PER cochannel receiver.
% Solve max sum log(1+a_i*p_i)-nu*sum(p_i) with 0<=p_i<=cap_i
% and sum(p_i)<=Pmax. nu>=0 is the shared C-UAV multiplier.
a=max(a(:),realmin);
cap=max(cap(:),0);
L=numel(a);
if L==0 || Pmax<=0
    p=zeros(L,1); return;
end
cap=min(cap,Pmax);
if sum(cap)<=Pmax
    if nu<=0
        p=cap;
    else
        p=min(cap,max(0,1/nu-1./a));
    end
    return;
end
if nu>0
    p0=min(cap,max(0,1/nu-1./a));
    if sum(p0)<=Pmax
        p=p0; return;
    end
end
lo=0; hi=1;
for t=1:120
    ph=min(cap,max(0,1/(nu+hi)-1./a));
    if sum(ph)<=Pmax, break; end
    hi=2*hi;
end
for t=1:120
    mid=(lo+hi)/2;
    pm=min(cap,max(0,1/(nu+mid)-1./a));
    if sum(pm)>Pmax, lo=mid; else, hi=mid; end
end
p=min(cap,max(0,1/(nu+hi)-1./a));
if sum(p)>Pmax, p=p*(Pmax/sum(p)); end
end
