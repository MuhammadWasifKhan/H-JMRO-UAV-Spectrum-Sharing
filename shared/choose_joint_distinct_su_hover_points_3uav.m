function [qHoverAll,bestIdx] = choose_joint_distinct_su_hover_points_3uav( ...
        qStartAll,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,area,seedBase,nAvg, ...
        useRayleigh,useShadow,useShadowA2G,Ntraj,stepMax)

% ICC ROBUST SU-ORIENTED TARGET RULE.

% L-UAV only: candidate targets must support a broad loop with half-span >=80 m.

% This affects initialization/target choice only; it is NOT a hard SCA shape constraint.

%

% The shortest flight time may not permit all three UAVs to visit three

% DISTINCT SU coordinates and return to their fixed start points.  Therefore:

%   1) Every UAV/SU pair first tries the exact SU position (alpha = 1).

%   2) If the exact SU is not reachable, the code finds the largest feasible

%      progress alpha in [0,1] along the straight direction from the UAV start

%      toward that SU.  This is a PHYSICAL flight-budget fallback, not an

%      arbitrary service-radius rule.

%   3) Three DISTINCT SU directions are then selected jointly.  The selection

%      first maximizes total reachable progress and then the role-aware

%      deterministic large-scale rate.

%

% Hence the code never aborts merely because T=200 s is too short for three

% exact SU visits.  As T increases, alpha naturally approaches 1 and exact SU

% service is used whenever feasible.

global UAV_HMIN UAV_HMAX UAV_DMIN_XY

U = 3;

K = size(WS,1);

altCand = unique([UAV_HMIN, qStartAll(1,3), UAV_HMAX]);

hsTry = [120 110 100 90 80 70 60 50 40 30 20];

feasible = false(U,K);

progress = zeros(U,K);       % alpha: 1 = exact SU, 0 = start point

score = -inf(U,K);

bestTarget = nan(U,K,3);

bestSpan = zeros(U,K);

startDist = inf(U,K);

for u = 1:U

    % Only the L-UAV was collapsing to a very narrow loop. For its target

    % feasibility test, reserve enough flight budget for a visibly broad

    % closed-loop initialization. C-UAV and U-UAV keep the original span set.

    if u == 2

        hsTryU = [120 110 100 90 80];

    else

        hsTryU = hsTry;

    end

    for k = 1:K

        q0 = qStartAll(u,:);

        suXY = WS(k,1:2);

        startDist(u,k) = norm(q0(1:2)-suXY);

        % ---------------------------------------------------------------

        % Find the largest reachable alpha toward this SU.

        % ---------------------------------------------------------------

        alphaBest = -1;

        zBest = NaN;

        spanBest = 0;

        % Exact target first.

        [ok,zTmp,spanTmp] = local_target_reachable( ...
            q0,suXY,1,altCand,Ntraj,stepMax,area,u,hsTryU);

        if ok

            alphaBest = 1;

            zBest = zTmp;

            spanBest = spanTmp;

        else

            % Bisection for the farthest feasible point toward the SU.

            % alpha=0 is the fixed start horizontal position and is always a

            % valid closed-loop fallback; a zero-width loop is permitted.

            lo = 0; hi = 1;

            [ok0,z0,sp0] = local_target_reachable( ...
                q0,suXY,0,altCand,Ntraj,stepMax,area,u,hsTryU);

            if ok0

                alphaBest = 0; zBest = z0; spanBest = sp0;

                for itB = 1:22

                    mid = 0.5*(lo+hi);

                    [okm,zm,spm] = local_target_reachable( ...
                        q0,suXY,mid,altCand,Ntraj,stepMax,area,u,hsTryU);

                    if okm

                        lo = mid;

                        alphaBest = mid;

                        zBest = zm;

                        spanBest = spm;

                    else

                        hi = mid;

                    end

                end

            end

        end

        if alphaBest < 0

            continue;

        end

        qXY = q0(1:2) + alphaBest*(suXY-q0(1:2));

        qc = [qXY,zBest];

        % Role-aware deterministic large-scale score at the physically

        % reachable target point.

        f = 0;

        for ss = 1:max(nAvg,1)

            rng(seedBase + 100000*u + 1000*k + ss,'twister');

            [huL,huU,hp2,hw2,hup,huw] = gen_channels_one_uav( ...
                qc,u,WS,WP,WW,w_b,w_wifi,H, ...
                f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
                useRayleigh,useShadow,useShadowA2G);

            if u == 1

                [rL,~,~,~] = allocate_one_band( ...
                    huL,hp2,P_PBS_j,sigma2_lic,hup,GammaLic,Pmax_lic,Bsc_lic);

                [rU,~,~,~] = allocate_one_band( ...
                    huU,hw2,P_WIFI_m,sigma2_un,huw,GammaUn,Pmax_un,Bsc_un);

                f = f + rL + rU;

            elseif u == 2

                [rL,~,~,~] = allocate_one_band( ...
                    huL,hp2,P_PBS_j,sigma2_lic,hup,GammaLic,Pmax_lic,Bsc_lic);

                f = f + rL;

            else

                [rU,~,~,~] = allocate_one_band( ...
                    huU,hw2,P_WIFI_m,sigma2_un,huw,GammaUn,Pmax_un,Bsc_un);

                f = f + rU;

            end

        end

        feasible(u,k) = true;

        progress(u,k) = alphaBest;

        score(u,k) = f/max(nAvg,1);

        bestTarget(u,k,:) = qc;

        bestSpan(u,k) = spanBest;

    end

end

exactCount = sum(progress >= 1-1e-5,2);

fprintf('  Exact BROAD-LOOP reachable SU counts at N=%d: C=%d, L=%d, U=%d\n', ...
    Ntraj,exactCount(1),exactCount(2),exactCount(3));

% -----------------------------------------------------------------------

% Jointly choose THREE DISTINCT SU directions.

% Primary: maximize total reachable progress alpha.

% Secondary: maximize deterministic role-aware rate.

% Tertiary: prefer smaller total start-to-SU distance.

% -----------------------------------------------------------------------

sepPrefs = [max(2*UAV_DMIN_XY,120), UAV_DMIN_XY, 0];

bestIdx = [NaN NaN NaN];

found = false;

for sepReq = sepPrefs

    bestProg = -inf;

    bestRate = -inf;

    bestDist = inf;

    for k1 = 1:K

        if ~feasible(1,k1), continue; end

        for k2 = 1:K

            if k2==k1 || ~feasible(2,k2), continue; end

            for k3 = 1:K

                if k3==k1 || k3==k2 || ~feasible(3,k3), continue; end

                pts = [squeeze(bestTarget(1,k1,1:2)).'; ...
                       squeeze(bestTarget(2,k2,1:2)).'; ...
                       squeeze(bestTarget(3,k3,1:2)).'];

                dTargets = [norm(pts(1,:)-pts(2,:)), ...
                            norm(pts(1,:)-pts(3,:)), ...
                            norm(pts(2,:)-pts(3,:))];

                if min(dTargets) < sepReq

                    continue;

                end

                pTot = progress(1,k1)+progress(2,k2)+progress(3,k3);

                rTot = score(1,k1)+score(2,k2)+score(3,k3);

                dTot = startDist(1,k1)+startDist(2,k2)+startDist(3,k3);

                if pTot > bestProg + 1e-9 || ...
                   (abs(pTot-bestProg)<=1e-9 && rTot > bestRate + 1e-6) || ...
                   (abs(pTot-bestProg)<=1e-9 && abs(rTot-bestRate)<=1e-6 && dTot < bestDist)

                    bestProg = pTot;

                    bestRate = rTot;

                    bestDist = dTot;

                    bestIdx = [k1 k2 k3];

                    found = true;

                end

            end

        end

    end

    if found, break; end

end

if ~found

    error(['Unable to assign three distinct SU directions. This should occur ' ...
           'only if K<3 or the node realization is invalid.']);

end

qHoverAll = zeros(U,3);

for u=1:U

    k=bestIdx(u);

    qHoverAll(u,:) = squeeze(bestTarget(u,k,:)).';

end

fprintf(['  Selected SU directions: C->SU%d (alpha=%.3f, span %.0f), ' ...
         'L->SU%d (alpha=%.3f, span %.0f), U->SU%d (alpha=%.3f, span %.0f)\n'], ...
    bestIdx(1),progress(1,bestIdx(1)),bestSpan(1,bestIdx(1)), ...
    bestIdx(2),progress(2,bestIdx(2)),bestSpan(2,bestIdx(2)), ...
    bestIdx(3),progress(3,bestIdx(3)),bestSpan(3,bestIdx(3)));

end
