function [rateL,rateU,idxL,pL,idxU,pU] = allocate_two_bands_shared_cuav( ...
        huL,hIL,PIL,sigL,hProtL,GammaL,PmaxL,BscL, ...
        huU,hIU,PIU,sigU,hProtU,GammaU,PmaxU,BscU,Ptotal)

% Joint one-to-one Hungarian assignment <-> exact KKT refinement for Mode C/C.

% Constraints:

%   sum(pL)<=PmaxL, sum(pU)<=PmaxU,

%   gL(j)*pL(j)<=GammaL for every j, gU(m)*pU(m)<=GammaU for every m,

%   sum(pL)+sum(pU)<=Ptotal.

global RESOURCE_AO_ITERS

if isempty(RESOURCE_AO_ITERS), RESOURCE_AO_ITERS=4; end

[K,J]=size(huL); [K2,M]=size(huU);

assert(K2==K,'Licensed/unlicensed user dimensions mismatch.');

pSubL=min(PmaxL/max(J,1),GammaL./max(hProtL(:),realmin));

pSubU=min(PmaxU/max(M,1),GammaU./max(hProtU(:),realmin));

s0=sum(pSubL)+sum(pSubU);

if s0>Ptotal && s0>0

    pSubL=pSubL*(Ptotal/s0); pSubU=pSubU*(Ptotal/s0);

end

prevPairsL=[]; prevPairsU=[];

idxL=zeros(0,2); idxU=zeros(0,2); pL=[]; pU=[];

for ao=1:max(1,RESOURCE_AO_ITERS)

    RmatL=zeros(K,J); RmatU=zeros(K,M);

    for k=1:K

        for j=1:J

            den=max(hIL(k,j)*PIL(j)+sigL,1e-30);

            RmatL(k,j)=BscL*log2(1+max(huL(k,j)*pSubL(j)/den,0));

        end

        for m=1:M

            den=max(hIU(k,m)*PIU(m)+sigU,1e-30);

            RmatU(k,m)=BscU*log2(1+max(huU(k,m)*pSubU(m)/den,0));

        end

    end

    pairsL=hungarian_one_to_one_max(RmatL);

    pairsU=hungarian_one_to_one_max(RmatU);

    idxL=pairsL; idxU=pairsU;

    [pL,pU]=exact_cuav_shared_power( ...
        idxL,huL,hIL,PIL,sigL,hProtL,GammaL,PmaxL, ...
        idxU,huU,hIU,PIU,sigU,hProtU,GammaU,PmaxU,Ptotal);

    pNewL=zeros(J,1); pNewU=zeros(M,1);

    for ii=1:size(idxL,1), pNewL(idxL(ii,2))=pL(ii); end

    for ii=1:size(idxU,1), pNewU(idxU(ii,2))=pU(ii); end

    samePairs=~isempty(prevPairsL) && isequal(prevPairsL,pairsL) && ...
              isequal(prevPairsU,pairsU);

    samePower=norm(pNewL-pSubL)<=1e-8*max(1,norm(pSubL)) && ...
              norm(pNewU-pSubU)<=1e-8*max(1,norm(pSubU));

    pSubL=pNewL; pSubU=pNewU;

    if samePairs && samePower, break; end

    prevPairsL=pairsL; prevPairsU=pairsU;

end

rateL=band_rate_from_assignment(idxL,pL,huL,hIL,PIL,sigL,BscL);

rateU=band_rate_from_assignment(idxU,pU,huU,hIU,PIU,sigU,BscU);

end
