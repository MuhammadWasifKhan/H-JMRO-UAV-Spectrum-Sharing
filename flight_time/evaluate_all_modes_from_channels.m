function [modeRates,modeInfo] = evaluate_all_modes_from_channels( ...
        huL,huU,hp2,hw2,hup,huw, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un)

modeTxL = [1,2,1,2];

modeTxU = [1,3,3,1];

modeRates = -inf(1,4);

modeInfo  = repmat(struct('txL',0,'txU',0,'rateL',0,'rateU',0, ...
    'txPowerL',0,'txPowerU',0,'txPowerTotal',0,'aggInterfL',0,'aggInterfU',0, ...
    'numActivePhysicalUAVs',0,'idxLic',[],'pLic',[],'idxUn',[],'pUn',[]),1,4);

global C_UAV_TOTAL_BUDGET

for mm=1:4

    txL=modeTxL(mm); txU=modeTxU(mm);

    hUL=squeeze(huL(:,:,txL));

    hUU=squeeze(huU(:,:,txU));

    if txL==1 && txU==1

        % Mode C/C: BOTH bands share the same physical C-UAV RF budget.

        % One-to-one Hungarian assignments and the joint shared-budget KKT powers are

        % alternated together, so the shared cap is part of the resource AO

        % rather than a post-hoc clipping step.

        [rateL,rateU,idxLic,pLic,idxUn,pUn] = allocate_two_bands_shared_cuav( ...
            hUL,hp2,P_PBS_j,sigma2_lic,hup(:,txL),GammaLic,Pmax_lic,Bsc_lic, ...
            hUU,hw2,P_WIFI_m,sigma2_un,huw(:,txU),GammaUn,Pmax_un,Bsc_un, ...
            C_UAV_TOTAL_BUDGET);

    else

        % Modes L/U, C/U and L/C: the two bands are carried by different

        % physical UAVs, hence each transmitting UAV uses its own band budget.

        [rateL,idxLic,pLic,~]=allocate_one_band( ...
            hUL,hp2,P_PBS_j,sigma2_lic,hup(:,txL),GammaLic,Pmax_lic,Bsc_lic);

        [rateU,idxUn,pUn,~]=allocate_one_band( ...
            hUU,hw2,P_WIFI_m,sigma2_un,huw(:,txU),GammaUn,Pmax_un,Bsc_un);

    end

    interfL=aggregate_interference(idxLic,pLic,hup(:,txL));

    interfU=aggregate_interference(idxUn,pUn,huw(:,txU));

    if any(cochannel_interference(idxLic,pLic,hup(:,txL))>GammaLic*(1+1e-7)) || ...
            any(cochannel_interference(idxUn,pUn,huw(:,txU))>GammaUn*(1+1e-7))

        error('Per-PU/WU cochannel interference constraint violated in mode %d.',mm);

    end

    if sum(pLic)>Pmax_lic*(1+1e-7) || sum(pUn)>Pmax_un*(1+1e-7)

        error('Band power constraint violated in mode %d.',mm);

    end

    if txL==1 && txU==1 && sum(pLic)+sum(pUn)>C_UAV_TOTAL_BUDGET*(1+1e-7)

        error('C-UAV shared total-power constraint violated.');

    end

    modeRates(mm)=rateL+rateU;

    modeInfo(mm).txL=txL; modeInfo(mm).txU=txU;

    modeInfo(mm).rateL=rateL; modeInfo(mm).rateU=rateU;

    modeInfo(mm).txPowerL=sum(pLic); modeInfo(mm).txPowerU=sum(pUn);

    modeInfo(mm).txPowerTotal=sum(pLic)+sum(pUn);

    modeInfo(mm).aggInterfL=interfL; modeInfo(mm).aggInterfU=interfU;

    modeInfo(mm).numActivePhysicalUAVs=numel(unique([txL,txU]));

    modeInfo(mm).idxLic=idxLic; modeInfo(mm).pLic=pLic;

    modeInfo(mm).idxUn=idxUn; modeInfo(mm).pUn=pUn;

end

end
