function [pL,pU] = exact_cuav_shared_power( ...
        idxL,huL,hIL,PIL,sigL,hProtL,GammaL,PmaxL, ...
        idxU,huU,hIU,PIU,sigU,hProtU,GammaU,PmaxU,Ptotal)

% Exact KKT solution of the fixed-assignment shared C-UAV power problem.

% No CVX call is used here, which is important because this function is called

% many times inside the Monte-Carlo/resource-AO loops.

LL=size(idxL,1); LU=size(idxU,1);

aL=zeros(LL,1); gL=zeros(LL,1); aU=zeros(LU,1); gU=zeros(LU,1);

for i=1:LL

    k=idxL(i,1); ss=idxL(i,2);

    aL(i)=max(huL(k,ss)/max(hIL(k,ss)*PIL(ss)+sigL,1e-30),1e-30);

    gL(i)=max(hProtL(ss),0);

end

for i=1:LU

    k=idxU(i,1); ss=idxU(i,2);

    aU(i)=max(huU(k,ss)/max(hIU(k,ss)*PIU(ss)+sigU,1e-30),1e-30);

    gU(i)=max(hProtU(ss),0);

end

% First solve with no shared-total multiplier. If the common cap is inactive,

% these are already the exact optimum.

pL=capped_waterfill(aL,GammaL./max(gL,realmin),PmaxL,0);

pU=capped_waterfill(aU,GammaU./max(gU,realmin),PmaxU,0);

if sum(pL)+sum(pU)<=Ptotal*(1+1e-10)

    return;

end

% Otherwise introduce the common KKT multiplier nu>=0. For a fixed nu the two

% bands separate; total power is monotone decreasing in nu, so bisection gives

% the exact common-budget KKT point.

nuLo=0; nuHi=1;

for kk=1:120

    pLt=capped_waterfill(aL,GammaL./max(gL,realmin),PmaxL,nuHi);

    pUt=capped_waterfill(aU,GammaU./max(gU,realmin),PmaxU,nuHi);

    if sum(pLt)+sum(pUt)<=Ptotal, break; end

    nuHi=2*nuHi;

end

for kk=1:120

    nuMid=0.5*(nuLo+nuHi);

    pLt=capped_waterfill(aL,GammaL./max(gL,realmin),PmaxL,nuMid);

    pUt=capped_waterfill(aU,GammaU./max(gU,realmin),PmaxU,nuMid);

    if sum(pLt)+sum(pUt)>Ptotal

        nuLo=nuMid;

    else

        nuHi=nuMid;

    end

end

pL=capped_waterfill(aL,GammaL./max(gL,realmin),PmaxL,nuHi);

pU=capped_waterfill(aU,GammaU./max(gU,realmin),PmaxU,nuHi);

% Tiny numerical clean-up only.

s=sum(pL)+sum(pU);

if s>Ptotal && s>0

    fac=Ptotal/s; pL=pL*fac; pU=pU*fac;

end

end
