function [p_assigned,idxList] = solve_band_kkt_aggregate( ...
        rho,hu2,hI,P_I_vec,sigma2,hProtect,Pmax,GammaAgg)

[~,S]=size(rho); idxList=zeros(0,2); a=[]; g=[];

for ss=1:S

    k=find(rho(:,ss)==1,1);

    if ~isempty(k)

        denom=max(hI(k,ss)*P_I_vec(ss)+sigma2,1e-30);

        idxList=[idxList;k,ss]; a=[a;max(hu2(k,ss)/denom,1e-30)]; g=[g;max(hProtect(ss),0)]; %#ok<AGROW>

    end

end

if isempty(a), p_assigned=[]; return; end

p_assigned=individual_capped_waterfill(a,g,Pmax,GammaAgg,0);

end
