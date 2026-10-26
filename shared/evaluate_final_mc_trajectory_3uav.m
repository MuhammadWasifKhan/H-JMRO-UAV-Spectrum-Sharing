function [avgRateMbps,avgTxPowerW,avgActiveUAV,modeFrac,slotMeanRateMbps] = ...
        evaluate_final_mc_trajectory_3uav( ...
        Qtraj,slotSeedIds,MC,baseSeed,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un, ...
        useRayleigh,useShadow,useShadowA2G)

% Final MC evaluator. It returns the MC-mean rate of each synchronized slot
% for diagnostics only. These values are NOT used to force monotonicity of
% the flight-time sweep and no fast-fading realization is duplicated across T.

[N,~,U] = size(Qtraj);

if numel(slotSeedIds) ~= N
    error('slotSeedIds length must equal the number of trajectory waypoints.');
end

totalRateMbps = 0;
totalTxPowerW = 0;
totalActive   = 0;
modeCount     = zeros(1,4);
slotRateAccum = zeros(1,N);

for mc = 1:MC
    sumRate = 0;
    sumPtx  = 0;
    sumAct  = 0;

    for n = 1:N
        rng(baseSeed + 1000000*mc + slotSeedIds(n),'twister');

        Qslot = zeros(U,3);
        for u = 1:U
            Qslot(u,:) = Qtraj(n,:,u);
        end

        [slotRate,bestMode,info] = eval_slot_3uav( ...
            Qslot,WS,WP,WW,w_b,w_wifi,H, ...
            f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
            P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
            GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
            Bsc_lic,Bsc_un,useRayleigh,useShadow,useShadowA2G);

        sumRate = sumRate + slotRate;
        sumPtx  = sumPtx + info.txPowerTotal;
        sumAct  = sumAct + info.numActivePhysicalUAVs;
        slotRateAccum(n) = slotRateAccum(n) + slotRate;
        modeCount(bestMode) = modeCount(bestMode) + 1;
    end

    totalRateMbps = totalRateMbps + (sumRate/N)/1e6;
    totalTxPowerW = totalTxPowerW + (sumPtx/N);
    totalActive   = totalActive + (sumAct/N);
end

avgRateMbps = totalRateMbps/MC;
avgTxPowerW = totalTxPowerW/MC;
avgActiveUAV = totalActive/MC;
modeFrac = modeCount/max(sum(modeCount),1);
slotMeanRateMbps = (slotRateAccum/MC)/1e6;

end
