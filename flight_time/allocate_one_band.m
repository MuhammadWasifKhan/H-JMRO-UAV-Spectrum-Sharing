function [sumRate,idxList,p_assigned,Rk] = allocate_one_band( ...
        hu2,hI,P_I_vec,sigma2,hProtect,GammaReceiver,Pmax,Bsc)

[K,S]=size(hu2);

assert(all(size(hI)==[K,S]),'hI must be K x S.');

assert(numel(P_I_vec)==S && numel(hProtect)==S,'Subchannel dimensions mismatch.');

global RESOURCE_AO_ITERS

if isempty(RESOURCE_AO_ITERS), RESOURCE_AO_ITERS=4; end

pSub=min(Pmax/max(S,1),GammaReceiver./max(hProtect(:),realmin)); prevPairs=[];

for ao=1:max(1,RESOURCE_AO_ITERS)

    Rmat=zeros(K,S);

    for k=1:K

        for ss=1:S

            denom=max(hI(k,ss)*P_I_vec(ss)+sigma2,1e-30);

            Rmat(k,ss)=Bsc*log2(1+max(hu2(k,ss)*max(pSub(ss),0)/denom,0));

        end

    end

    pairs=hungarian_one_to_one_max(Rmat); rho=zeros(K,S);

    for ii=1:size(pairs,1), rho(pairs(ii,1),pairs(ii,2))=1; end

    [p_assigned,idxList]=solve_band_kkt_per_receiver(rho,hu2,hI,P_I_vec,sigma2,hProtect,Pmax,GammaReceiver);

    pNew=zeros(S,1); for ell=1:size(idxList,1), pNew(idxList(ell,2))=p_assigned(ell); end

    if ~isempty(prevPairs) && isequal(prevPairs,pairs) && norm(pNew-pSub)<=1e-8*max(1,norm(pSub)), pSub=pNew; break; end

    prevPairs=pairs; pSub=pNew;

end

Rk=zeros(K,1); sumRate=0;

for ell=1:size(idxList,1)

    k=idxList(ell,1); ss=idxList(ell,2); denom=max(hI(k,ss)*P_I_vec(ss)+sigma2,1e-30);

    rr=Bsc*log2(1+max(hu2(k,ss)*p_assigned(ell)/denom,0)); Rk(k)=Rk(k)+rr; sumRate=sumRate+rr;

end

if sum(p_assigned)>Pmax*(1+1e-8), error('KKT total-power violation.'); end

if any(cochannel_interference(idxList,p_assigned,hProtect)>GammaReceiver*(1+1e-8)), error('KKT per-receiver interference violation.'); end

end
