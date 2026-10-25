function rate = band_rate_from_assignment(idxList,p,hu2,hI,P_I_vec,sigma2,Bsc)

rate=0;

for ell=1:size(idxList,1)

    k=idxList(ell,1); ss=idxList(ell,2); denom=max(hI(k,ss)*P_I_vec(ss)+sigma2,1e-30);

    rate=rate+Bsc*log2(1+max(hu2(k,ss)*max(p(ell),0)/denom,0));

end

end
