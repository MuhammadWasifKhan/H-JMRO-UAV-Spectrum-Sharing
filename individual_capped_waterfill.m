function p = individual_capped_waterfill(a,g,Pmax,Gamma,nu)
% Max sum_i log(1+a_i*p_i)-nu*sum_i p_i
% s.t. sum_i p_i<=Pmax and g_i*p_i<=Gamma, independently for each i.
% For fixed nu, the per-band total-power multiplier is found by bisection.
a = max(a(:),1e-30);
g = max(g(:),0);
L = numel(a);
p = zeros(L,1);
if L == 0 || Pmax <= 0 || Gamma <= 0, return; end
cap = Gamma ./ max(g,realmin);
cap(g==0) = Pmax;
cap = min(cap,Pmax);
if nu > 0
    p0 = min(cap,max(0,1/nu-1./a));
    if sum(p0) <= Pmax*(1+1e-12), p=p0; return; end
end
lo = 0; hi = 1;
for j=1:120
    trial = min(cap,max(0,1/(hi+nu)-1./a));
    if sum(trial) <= Pmax, break; end
    hi = hi*2;
end
for j=1:100
    mid = (lo+hi)/2;
    trial = min(cap,max(0,1/(mid+nu)-1./a));
    if sum(trial) > Pmax, lo=mid; else, hi=mid; end
end
p = min(cap,max(0,1/(hi+nu)-1./a));
end
