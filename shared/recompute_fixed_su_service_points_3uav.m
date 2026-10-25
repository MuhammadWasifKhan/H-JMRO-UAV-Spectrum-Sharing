function qHoverAll = recompute_fixed_su_service_points_3uav( ...
        targetSUIdx,qStartAll,WS,Ntraj,stepMax,area)

% Keep SU identities fixed across the T sweep. For each current flight time,

% find the farthest broad-loop-feasible service point along start -> fixed SU.

% This changes only reachability with T; it never switches the assigned SU.

U = 3;

qHoverAll = zeros(U,3);

for u=1:U

    q0 = qStartAll(u,:);

    k  = targetSUIdx(u);

    qSU = [WS(k,1), WS(k,2), q0(3)];

    qSU = clip_q3d_role(qSU,area,u);

    % Exact SU first.

    okExact = false;

    for hs = [120 110 100 90 80 70 60 50 40 30 20]

        [~,ok] = make_broad_service_loop(q0,qSU,Ntraj,stepMax,area,u,hs);

        if ok

            qHoverAll(u,:) = qSU;

            okExact = true;

            break;

        end

    end

    if okExact

        continue;

    end

    % Otherwise use the farthest feasible point toward the SAME SU.

    lo = 0; hi = 1; qBest = q0;

    for it=1:50

        a = 0.5*(lo+hi);

        qTry = clip_q3d_role(q0 + a*(qSU-q0),area,u);

        okAny = false;

        for hs = [120 110 100 90 80 70 60 50 40 30 20]

            [~,ok] = make_broad_service_loop(q0,qTry,Ntraj,stepMax,area,u,hs);

            if ok

                okAny = true;

                break;

            end

        end

        if okAny

            lo = a;

            qBest = qTry;

        else

            hi = a;

        end

    end

    qHoverAll(u,:) = qBest;

end

end

function [avgRateMbps,avgTxPowerW,avgActiveUAV,modeFrac,slotMeanRateMbps] = ...
        evaluate_final_mc_trajectory_3uav( ...
        Qtraj,slotSeedIds,MC,baseSeed,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un, ...
        useRayleigh,useShadow,useShadowA2G)

% Final MC evaluator.  It returns the MC-mean rate of each synchronized slot

% for diagnostics only.  These values are NOT used to force monotonicity of

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

        for u=1:U

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
