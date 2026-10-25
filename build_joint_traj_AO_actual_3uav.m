function [Qout,convHistory] = build_joint_traj_AO_actual_3uav( ...
        Qinit,q_start_all,q_end_all,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,AO_traj_iter,stepMax,alpha_traj,area, ...
        seedBase,trajEvalSamples,useRayleigh,useShadow,useShadowA2G)

% =========================================================================
% JOURNAL-STYLE AO-SCA TRAJECTORY FOR THE HETEROGENEOUS C/L/U ICC SCHEME
% =========================================================================
% This uses the SAME trajectory logic as the working 2-C-UAV comparison:
%   1) Re-optimize mode selection + Hungarian allocation + KKT power at Q^(r).
%   2) Freeze those communication variables only to form the analytical
%      first-order rate gradient of the Al-Hourani channel.
%   3) Solve ONE joint CVX trajectory subproblem over all three UAVs and all
%      waypoints using geometry, mobility, endpoint, separation and trust-region
%      constraints. Frozen-power PU/WU Taylor constraints are NOT duplicated
%      inside CVX.
%   4) At every candidate trajectory, re-solve the COMPLETE communication
%      problem. Therefore the original aggregate PU/WU interference protection
%      is re-enforced by the normal mode/Hungarian/KKT resource optimizer.
%   5) Accept the candidate only when the actual deterministic re-optimized
%      average sum rate is nondecreasing.
%
% The heterogeneous architecture itself is unchanged: C/C, L/U, C/U and L/C
% mode selection, one transmitter per band, Hungarian assignment and KKT power.
% =========================================================================

global UAV_HMIN UAV_HMAX UAV_DMIN UAV_DMIN_XY UAV_XMIN UAV_XMAX

[N,~,U] = size(Qinit);
Qout = Qinit;
trajUseRayleigh = false;

% Make the initial trajectory exactly feasible before starting SCA.
for u = 1:U
    Qout(:,:,u) = enforce_trajectory_feasibility_single( ...
        Qout(:,:,u),area,stepMax,q_start_all(u,:),q_end_all(u,:),u);
end

Qout = repair_joint_separation_3uav( ...
    Qout,area,stepMax,q_start_all,q_end_all,UAV_DMIN);

if min_continuous_segment_separation_3uav(Qout) < UAV_DMIN*(1-1e-8) || ...
   min_continuous_segment_horizontal_separation_3uav(Qout) < UAV_DMIN_XY*(1-1e-8)
    warning(['Initial trajectories are not continuously separation-feasible. ' ...
             'Using stationary region-feasible closed-loop initialization.']);
    for u = 1:U
        Qout(:,:,u) = repmat(q_start_all(u,:),N,1);
    end
end

if min_continuous_segment_separation_3uav(Qout) < UAV_DMIN*(1-1e-8) || ...
   min_continuous_segment_horizontal_separation_3uav(Qout) < UAV_DMIN_XY*(1-1e-8)
    error('The UAV deployment/start points violate the required separation.');
end

% Exact reference-feasibility audit BEFORE SCA.  Never start an SCA iteration
% from a point that violates the same geometry tests used for acceptance.
for u = 1:U
    Qout(:,:,u) = enforce_trajectory_feasibility_single( ...
        Qout(:,:,u),area,stepMax,q_start_all(u,:),q_end_all(u,:),u);
    st0 = sqrt(sum(diff(Qout(:,:,u),1,1).^2,2));
    if max(st0) > stepMax*(1+1e-8) || ...
       any(Qout(:,1,u) < UAV_XMIN(u)-1e-9) || ...
       any(Qout(:,1,u) > UAV_XMAX(u)+1e-9) || ...
       any(Qout(:,2,u) < -1e-9) || any(Qout(:,2,u) > area+1e-9) || ...
       any(Qout(:,3,u) < UAV_HMIN-1e-9) || ...
       any(Qout(:,3,u) > UAV_HMAX+1e-9)
        error('Initial/reference trajectory is not geometry-feasible for UAV %d.',u);
    end
end

% No fixed reference-loop corridor: the SCA trust region controls each local step.
if exist('cvx_begin','file') ~= 2
    error('CVX is required to compute actual AO-SCA convergence.');
end

trust0 = max(1e-3,min(stepMax,alpha_traj*stepMax));
trustMin = max(0.05,trust0/64);
maxRetries = 7;
scaTol = 1e-4;

% Deterministic large-scale-CSI objective; the COMPLETE communication problem
% is solved by evaluate_trajectory_actual_3uav at every accepted iterate.
fPrev = evaluate_trajectory_actual_3uav( ...
    Qout,WS,WP,WW,w_b,w_wifi,H, ...
    f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
    P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
    GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
    Bsc_lic,Bsc_un,seedBase+999,trajEvalSamples, ...
    trajUseRayleigh,useShadow,useShadowA2G);

convHistory = nan(1,AO_traj_iter+1);
convHistory(1) = fPrev/1e6;

for ell = 1:AO_traj_iter

    Qbase = Qout;

    % Re-solve mode + Hungarian + KKT at the CURRENT trajectory and construct
    % the exact analytical first-order frozen-resource rate gradient.
    scaData = cell(N,1);
    fixedInfo = cell(N,1);
    for n = 1:N
        Qslot = squeeze(Qbase(n,:,:)).';
        [scaData{n},fixedInfo{n}] = build_reference_sca_data_3uav( ...
            Qslot,WS,WP,WW,w_b,w_wifi,H, ...
            f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
            P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
            GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
            Bsc_lic,Bsc_un,useShadow,useShadowA2G);
    end

    accepted = false;
    acceptedRate = fPrev;
    radius = trust0;

    for retry = 1:maxRetries

        if radius < trustMin
            break;
        end

        try
            [Qcand,statusNow] = solve_traj_cvx_subproblem_3uav( ...
                Qbase,scaData,q_start_all,q_end_all,N,U, ...
                radius,stepMax,0,GammaLic,GammaUn,area);
        catch ME
            clear cvx_problem
            fprintf(['  SCA iter %d/%d | retry %d/%d | trust=%.4f m | ' ...
                     'CVX error: %s\n'], ...
                ell,AO_traj_iter,retry,maxRetries,radius,ME.message);
            if exist('cvx_clear','file') == 2
                try
                    cvx_clear;
                catch
                end
            end
            radius = radius/2;
            continue;
        end

        fprintf('  SCA iter %d/%d | retry %d/%d | trust=%.4f m | CVX=%s\n', ...
            ell,AO_traj_iter,retry,maxRetries,radius,statusNow);

        solved = contains(lower(statusNow),'solved') || ...
                 contains(lower(statusNow),'inaccurate');

        if solved
            % CVX solutions can have small numerical residuals in SOC mobility
            % constraints.  Do not accept such a point blindly.  Backtrack from
            % the CVX candidate toward the already-feasible expansion point and
            % accept only an EXACTLY geometry-feasible candidate.
            for bt = 0:12
                etaStep = 0.5^bt;
                Qt = Qbase + etaStep*(Qcand-Qbase);

                % Repair the backtracked candidate BEFORE any acceptance test.
                % This removes CVX numerical residues in box/altitude/mobility
                % constraints while preserving the closed-loop endpoints.
                % IMPORTANT: the repaired point itself is subsequently checked
                % for mobility, separation and objective; no
                % post-acceptance clipping/projection is performed.
                for uu = 1:U
                    Qt(:,:,uu) = enforce_trajectory_feasibility_single( ...
                        Qt(:,:,uu),area,stepMax,q_start_all(uu,:), ...
                        q_end_all(uu,:),uu);
                end

                geomOK = true;
                for uu = 1:U
                    st = sqrt(sum(diff(Qt(:,:,uu),1,1).^2,2));
                    if max(st) > stepMax*(1+1e-6)
                        geomOK = false;
                        break;
                    end
                    if any(Qt(:,1,uu) < UAV_XMIN(uu)-1e-9) || ...
                       any(Qt(:,1,uu) > UAV_XMAX(uu)+1e-9) || ...
                       any(Qt(:,2,uu) < -1e-9) || ...
                       any(Qt(:,2,uu) > area+1e-9) || ...
                       any(Qt(:,3,uu) < UAV_HMIN-1e-9) || ...
                       any(Qt(:,3,uu) > UAV_HMAX+1e-9)
                        geomOK = false;
                        break;
                    end
                end

                if ~geomOK
                    continue;
                end

                sep3D = min_continuous_segment_separation_3uav(Qt);
                sepXY = min_continuous_segment_horizontal_separation_3uav(Qt);
                % Keep a 1 mm feasibility margin against CVX floating-point
                % residuals. A candidate must genuinely satisfy Dmin=100 m.
                % The last accepted trajectory remains unchanged if no
                % sufficiently safe improving candidate is found.
                sepSafetyMargin = 1e-3;  % metres
                if sep3D < UAV_DMIN + sepSafetyMargin || ...
                   sepXY < UAV_DMIN_XY + sepSafetyMargin
                    continue;
                end

                % Validate original nonlinear per-PU/WU interference with
                % frozen reference powers before reoptimizing resources.
                if ~fixed_resource_interference_feasible_3uav( ...
                        Qt,fixedInfo,WP,WW,f_lic,f_un,c,A2G,Shadow, ...
                        GammaLic,GammaUn,useShadowA2G)
                    continue;
                end

                % Re-optimize the COMPLETE communication problem at the exact
                % geometry-feasible candidate, as in the working comparison logic.
                fCand = evaluate_trajectory_actual_3uav( ...
                    Qt,WS,WP,WW,w_b,w_wifi,H, ...
                    f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
                    P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
                    GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
                    Bsc_lic,Bsc_un,seedBase+999,trajEvalSamples, ...
                    trajUseRayleigh,useShadow,useShadowA2G);

                if isfinite(fCand) && ...
                   fCand >= fPrev - 1e-9*max(1,abs(fPrev))
                    relImprove = abs(fCand-fPrev)/max(1,abs(fPrev));
                    Qout = Qt;
                    fPrev = fCand;
                    acceptedRate = fCand;
                    accepted = true;

                    fprintf(['    accepted | step=%.4f | reoptimized deterministic ' ...
                             'rate = %.6f Mbps | relative gain = %.3e\n'], ...
                        etaStep,acceptedRate/1e6,relImprove);
                    break;
                end
            end

            if accepted
                break;
            end
        end

        radius = radius/2;
    end

    convHistory(ell+1) = fPrev/1e6;

    if ~accepted
        fprintf(['  No feasible nondecreasing journal-style SCA step; ' ...
                 'trajectory iteration stops.\n']);
        break;
    end

    if relImprove <= scaTol
        break;
    end
end

% Carry the last *actual* accepted objective forward after termination.
% The resulting plateau reflects algorithm termination, not fabricated gains.
lastRecorded = find(isfinite(convHistory),1,'last');
convHistory(lastRecorded+1:end) = convHistory(lastRecorded);

% Final exact geometry checks.
for u = 1:U
    Qout(1,:,u) = q_start_all(u,:);
    Qout(end,:,u) = q_end_all(u,:);

    steps = sqrt(sum(diff(Qout(:,:,u),1,1).^2,2));
    if max(steps) > stepMax*(1+1e-6)
        error('Final SCA trajectory violates speed constraint for UAV %d.',u);
    end

    if any(Qout(:,1,u) < UAV_XMIN(u)-1e-8) || ...
       any(Qout(:,1,u) > UAV_XMAX(u)+1e-8)
        error('Final SCA trajectory violates the horizontal flight area for UAV %d.',u);
    end

    if any(Qout(:,2,u) < -1e-8) || any(Qout(:,2,u) > area+1e-8) || ...
       any(Qout(:,3,u) < UAV_HMIN-1e-8) || any(Qout(:,3,u) > UAV_HMAX+1e-8)
        error('Final SCA trajectory violates y/altitude bounds for UAV %d.',u);
    end
end

if min_joint_separation_3uav(Qout) < UAV_DMIN*(1-1e-8)
    error('Final SCA trajectory violates waypoint inter-UAV separation.');
end

if min_continuous_segment_separation_3uav(Qout) < UAV_DMIN*(1-1e-8)
    error(['Final SCA trajectory violates continuous 3D inter-UAV separation ' ...
           'between consecutive waypoints.']);
end

if min_continuous_segment_horizontal_separation_3uav(Qout) < UAV_DMIN_XY*(1-1e-8)
    error(['Final SCA trajectory violates continuous horizontal inter-UAV ' ...
           'separation between consecutive waypoints.']);
end
end
