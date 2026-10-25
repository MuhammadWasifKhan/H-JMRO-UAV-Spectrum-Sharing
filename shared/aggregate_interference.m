function I = aggregate_interference(idxList,p,hProtect)

I = 0;

for ell=1:size(idxList,1)

    ss = idxList(ell,2);

    I = I + max(p(ell),0)*max(hProtect(ss),0);

end

end
