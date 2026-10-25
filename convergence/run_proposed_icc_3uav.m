function proposedRates = run_proposed_icc_3uav()
% =========================================================================

% ICC PROPOSED SCHEME -- SU-ORIENTED CLOSED-LOOP / NON-COLLIDING VERSION

% INTERFERENCE-AWARE HETEROGENEOUS 3-UAV MODE SELECTION + EXACT-TAYLOR 3D TRAJECTORY

%

% UAV-1 : Cognitive UAV (C-UAV)       -> licensed + unlicensed capable

% UAV-2 : Licensed-only UAV (L-UAV)   -> licensed only

% UAV-3 : Unlicensed-only UAV (U-UAV) -> unlicensed only

%

% PROPOSED H-JMRO ARCHITECTURE:

%   1) Unlicensed access is NOT QoS-triggered.

%   2) Licensed and unlicensed bands are considered in EVERY slot.

%   3) At each synchronized slot, the controller tests four feasible modes:

%

%      Mode 1: C-UAV on licensed  + C-UAV on unlicensed   = C/C

%      Mode 2: L-UAV on licensed  + U-UAV on unlicensed   = L/U

%      Mode 3: C-UAV on licensed  + U-UAV on unlicensed   = C/U

%      Mode 4: L-UAV on licensed  + C-UAV on unlicensed   = L/C

%

%   4) For every candidate mode:

%         iterative one-to-one Hungarian assignment <-> per-receiver-interference KKT power refinement

%            -> post-power sum-rate calculation

%   5) MODE CONTROL uses instantaneous composite CSI: G = beta*|g|^2.

%      In each slot, one band-specific channel realization is generated for every

%      physical UAV-receiver link and reused consistently across all four candidate

%      modes. Different UAVs and different bands have different channel realizations.

%   6) All three UAV trajectories are synchronized and optimized in 3D by

%      exact first-order Taylor SCA/SCP of the Al-Hourani model. All three UAVs

%      can move across the full service area; no role-specific corridor is imposed.

%      At each iteration

%      the selected mode, subchannel assignment, and powers are fixed.

%   7) Desired-rate and individual interference Taylor gradients are obtained

%      analytically from the full Al-Hourani A2G model; no finite differences are used.

%      The resulting affine local surrogate is solved jointly over all waypoints

%      and all three UAVs using CVX under mobility, endpoint and trust constraints.

%   8) Each SCA candidate is checked against the original nonlinear fixed-resource

%      interference constraints and actual deterministic objective. Adaptive

%      trust-region re-solving plus backtracking preserves feasibility and

%      nondecreasing accepted outer iterates.

%   9) UAV-to-SU and ground-interferer-to-SU Rayleigh fading is explicitly

%      subchannel dependent. Instantaneous composite CSI is used for mode/resource

%      selection, while trajectory SCA remains based on deterministic large-scale

%      Al-Hourani CSI. A minimum 3-D inter-UAV separation is enforced by SCA.

%  10) All three UAVs may use the complete 1000 m x 1000 m horizontal area.

%      Collision avoidance is time synchronized: 3-D and horizontal separation

%      are enforced at waypoints and continuously between consecutive slots.

%      Geometric path crossings at DIFFERENT times are allowed.

%  11) Each mission is initialized as a genuine non-retracing closed loop:

%      start -> broad curved branch -> SU-oriented service point -> broad return branch -> start.

%      The exact endpoint equality is also enforced inside every SCA subproblem.

%  12) QoS is intentionally not imposed: the objective is sum-rate maximization

%      under power, mobility, separation, and individual cochannel interference protection.

%      Mode-switching latency/coordination overhead is assumed negligible relative

%      to the slot duration and is therefore not charged in the objective.

%

% Thus there is ONE transmitter per band in every slot. Depending on the

% selected mode, either one physical UAV (Mode 1) or two physical UAVs

% (Modes 2-4) transmit in that slot.

%

% NOTE:

%   "Communication-only EE" reported below includes RF transmit power only.

%   Propulsion/hovering power is intentionally NOT included in that metric.

% =========================================================================

if exist('cvx_clear','file') == 2

    try

        cvx_clear;

    catch

    end

end

%% ========================== System parameters ============================

area = 1000;

% TRUE 3D altitude settings for all three UAVs.

H0   = 100;      % initial/final altitude (m)

Hmin = 50;       % minimum altitude (m)

Hmax = 150;      % maximum altitude (m)

% H is retained only for legacy function signatures. The actual UAV altitude

% used in every A2G channel calculation is q(3).

H = H0;

global UAV_HMIN UAV_HMAX RESOURCE_AO_ITERS C_UAV_TOTAL_BUDGET UAV_DMIN UAV_DMIN_XY UAV_XMIN UAV_XMAX

UAV_HMIN = Hmin;

UAV_HMAX = Hmax;

K = 6;       % secondary users

J = 4;       % licensed subchannels / primary users

M = 4;       % unlicensed subchannels / Wi-Fi users

% B_total_lic = 10e6;

% B_total_un  = 10e6;

Bsc_lic = 180e3;

Bsc_un  = 180e3;

%% ========================== Noise powers ================================

N0_dBmHz = -169;

sigma2_lic_dBm = N0_dBmHz + 10*log10(Bsc_lic);

sigma2_un_dBm  = N0_dBmHz + 10*log10(Bsc_un);

sigma2_lic = 10.^((sigma2_lic_dBm-30)/10);

sigma2_un  = 10.^((sigma2_un_dBm -30)/10);

%% ========================== UAV power budgets ===========================

% Independent per-band budgets.

% For the C-UAV, these are its licensed and unlicensed budgets.

% For L-UAV / U-UAV, only the relevant band is used.

Pmax_lic = 0.5;  % W per licensed-band transmitter

Pmax_un  = 0.5;  % W per unlicensed-band transmitter

% Overall RF hardware cap for the cognitive UAV when it serves BOTH bands.

% With 0.5 W per-band limits and a 1.0 W total cap, this total constraint is

% non-binding by itself, but it is retained explicitly for hardware accounting.

Pmax_C_total = 1.0;  % W, C-UAV total RF cap in Mode C/C

%% ========================== External interferers ========================

P_PBS  = 2.5;    % W per licensed subchannel

P_WIFI = 1.0;    % W per unlicensed subchannel

P_PBS_j  = P_PBS  * ones(J,1);

P_WIFI_m = P_WIFI * ones(M,1);

%% ========================== Carrier frequencies =========================

c = 3e8;

f_lic  =2.0e9;

f_un   = 5.0e9;

f_pbs  = f_lic;

f_wifi = f_un;

%% ========================== Propagation =================================

n_pbs  = 3.7;

n_wifi = 3.0;

sigmaSh_pbs_dB  = 8;

sigmaSh_wifi_dB = 6;

A2G.a = 9.61;

A2G.b = 0.16;

A2G.etaLoS_dB    = 1.0;

A2G.etaNLoS_dB   = 20.0;

A2G.sigmaLoS_dB  = 4;

A2G.sigmaNLoS_dB = 6;

%% ========================== Interference limits =========================

GammaLic_dBm = -80;

GammaUn_dBm  = -80;

GammaLic = 10.^((GammaLic_dBm-30)/10);  % W, individual interference limit per licensed PU/subchannel

GammaUn  = 10.^((GammaUn_dBm -30)/10);  % W, individual interference limit per Wi-Fi user/subchannel

%% ========================== Monte Carlo =================================

MC     = 300;     % fast-fading realizations per geometry

NDrops = 30;     % independent spatial drops for statistically robust results

dropSeedBase = 9000;

%% ========================== Ground nodes ================================

w_b    = [area/2, area/2];

w_wifi = [area/2, area/2 + 200];

%% ========================== Three synchronized UAVs =====================

U = 3;

uavName = {'C-UAV','L-UAV','U-UAV'};

% All three heterogeneous UAVs may move over the complete horizontal area.

% Code order is: 1=C-UAV, 2=L-UAV, 3=U-UAV.

% This permits natural optimized loops and also permits two XY path traces to

% cross geometrically when the UAVs reach the crossing point at different times.

% Collision safety is enforced by synchronized continuous-time separation.

UAV_XMIN = [0, 0, 0];

UAV_XMAX = [area, area, area];

% Distinct start/end points; missions remain closed loops.

% q_start_all = [ ...
%     500, 500, H0;   % C-UAV
% 
%     150, 250, H0;   % L-UAV
% 
%     850, 750, H0];  % U-UAV500, 500, H0
%%%%%%%%%%%%120 850;260 180;430 620
q_start_all = [ ...
    500.0, 600.0, H0;   % C-UAV
    413.4, 450.0, H0;   % L-UAV
    586.6, 450.0, H0];  % U-UAV

q_end_all = q_start_all;

Dmin = 100;   % minimum inter-UAV separation (m)
% q_start_all = [ ...
%     500.00, 520.00, H0;   % C-UAV 1
%     482.68, 490.00, H0;   % C-UAV 2
%     517.32, 490.00, H0];  % C-UAV 3
% 
% % q_end_all = q_start_all;
% % q_start_all = [ ...
% %     200, 800, H0;   % C-UAV
% % 
% %     500, 500, H0;   % L-UAV
% % 
% %     800, 200, H0];  % U-UAV
% 
% q_end_all = q_start_all;

Vmax_user = 5;             % m/s

UAV_DMIN = Dmin;           % 100 m minimum synchronized continuous 3D separation

UAV_DMIN_XY = Dmin;        % 100 m minimum synchronized continuous XY separation

for u=1:U

    assert(q_start_all(u,1) >= 0 && q_start_all(u,1) <= area && ...
           q_start_all(u,2) >= 0 && q_start_all(u,2) <= area, ...
        '%s start point is outside the flight area.',uavName{u});

end

% Verify 100 m separation at the fixed start/end points before initialization.
for u = 1:U-1
    for v = u+1:U
        d3 = norm(q_start_all(u,:) - q_start_all(v,:));
        dxy = norm(q_start_all(u,1:2) - q_start_all(v,1:2));
        assert(d3 >= UAV_DMIN - 1e-8 && dxy >= UAV_DMIN_XY - 1e-8, ...
            'Initial positions of %s and %s violate 100 m separation.', ...
            uavName{u},uavName{v});
    end
end

%% ========================== Trajectory controls =========================

AO_traj_iter  = 40;

RESOURCE_AO_ITERS = 3;       % equivalent to AO_AP_iter = 3 in journal code

C_UAV_TOTAL_BUDGET = Pmax_C_total;

alpha_traj    = 0.30;   % local SCA trust-region scale

trajEvalSamples = 1;   % trajectory design uses deterministic large-scale channels

trajSeedBase  = 77777;

%% ========================== Randomness ==================================

useRayleigh  = true;

useShadow    = true;

useShadowA2G = true;

%% ========================== Flight-time sweep ===========================

T_list = 220;

dt_fixed = 5;

STEP_MAX = Vmax_user*dt_fixed;

assert(all(abs(T_list/dt_fixed-round(T_list/dt_fixed)) < 1e-12), ...
    'Every T must be an integer multiple of dt_fixed.');

fprintf('\n============================================================\n');

fprintf('ICC HETEROGENEOUS 3-UAV INTERFERENCE-AWARE SCHEME\n');

fprintf('============================================================\n');

fprintf('UAVs: C-UAV + L-UAV + U-UAV\n');

fprintf('Both bands considered in every synchronized slot; no QoS constraint.\n');

fprintf('Modes: C/C, L/U, C/U, L/C\n');

fprintf('Inter-UAV cochannel communication interference: NONE (one UAV transmitter per band in every mode)\n');

fprintf('T-sweep logic: fixed SU identities per drop + nested previous-T warm start + paired CRN (no forced monotonicity)\n');

fprintf('K=%d, J=%d, M=%d, altitude range=[%.0f, %.0f] m\n',K,J,M,Hmin,Hmax);

fprintf('Pmax_lic=%.2f W, Pmax_un=%.2f W, C-UAV total cap=%.2f W\n',Pmax_lic,Pmax_un,Pmax_C_total);

fprintf('GammaLic=%.0f dBm, GammaUn=%.0f dBm\n',GammaLic_dBm,GammaUn_dBm);

fprintf('dt=%.2f s, Vmax=%.2f m/s, STEP_MAX=%.2f m\n', ...
    dt_fixed,Vmax_user,STEP_MAX);

fprintf('Resource AO iterations per mode = %d\n',RESOURCE_AO_ITERS);

fprintf('Minimum inter-UAV separation = %.1f m (3D) and %.1f m synchronized XY\n',UAV_DMIN,UAV_DMIN_XY);

fprintf('Horizontal mobility: all UAVs use the full [0,1000] x [0,1000] m area\n');

fprintf('Geometric XY path crossings are allowed only when they occur at different times\n');

fprintf('Trajectory targeting: first-T reachable SU-oriented closed-loop initialization; free SCA thereafter\n');

fprintf('No UAV-specific trajectory corridors are imposed.\n');

fprintf('Trajectory method: exact first-order Taylor SCA/SCP of full Al-Hourani model\n');

fprintf('CVX installed? exist(cvx_begin) = %d\n\n',exist('cvx_begin','file'));

validate_a2g_gradient_once(A2G,c,f_lic,useShadowA2G);

%% ========================== Convergence history ========================
convTraceDrop = nan(NDrops,AO_traj_iter+1,numel(T_list));

%% ========================== Accumulators ================================

accRateDrop     = zeros(1,AO_traj_iter);

accTxPowerDrop  = zeros(size(T_list));

accActiveDrop   = zeros(size(T_list));

accModeFracDrop = zeros(numel(T_list),4);

%% ========================================================================

%                                DROP LOOP

% ========================================================================

for drop = 1:NDrops

    rng(dropSeedBase + drop,'twister');

    % Fixed node realization for a drop.

    WS = area*rand(K,2);

    WP = area*rand(J,2);

    WW = area*rand(M,2);

    % Coherent shadowing retained over the drop.

    Shadow = struct();

    if useShadow

        Shadow.pbs_su_dB  = sigmaSh_pbs_dB*randn(K,1);

        Shadow.wifi_su_dB = sigmaSh_wifi_dB*randn(K,1);

    else

        Shadow.pbs_su_dB  = zeros(K,1);

        Shadow.wifi_su_dB = zeros(K,1);

    end

    if useShadowA2G

        Shadow.a2g_su_z = randn(K,U);

        Shadow.a2g_pu_z = randn(J,U);

        Shadow.a2g_wu_z = randn(M,U);

    else

        Shadow.a2g_su_z = zeros(K,U);

        Shadow.a2g_pu_z = zeros(J,U);

        Shadow.a2g_wu_z = zeros(M,U);

    end

    baseSeed = 12345 + 100000*drop;

    avgRate     = zeros(size(T_list));

    avgTxPower  = zeros(size(T_list));

    avgActive   = zeros(size(T_list));

    modeFrac    = zeros(numel(T_list),4);

    Qprev = [];

    slotSeedPrev = [];   % RNG identities attached to inherited old slots (CRN only)

    %% ======================== Flight-time sweep =========================

    for t_idx = 1:numel(T_list)

        T_fixed = T_list(t_idx);

        N = round(T_fixed/dt_fixed) + 1;

        if abs((N-1)*dt_fixed - T_fixed) > 1e-9

            error('T and dt_fixed are inconsistent.');

        end

        % ---------------------------------------------------------------

        % 1) JOINT FLIGHT-BUDGET-CONSTRAINED SU-oriented service-point selection.

        % ---------------------------------------------------------------

        % IMPORTANT FIX:

        % The controller jointly selects THREE DISTINCT SU directions (one per

        % physical UAV). Exact SU visits are preferred. If the current T is too

        % short, the target is the farthest physically reachable point toward

        % that SU, so the simulation remains feasible without an arbitrary

        % service-radius rule.

        hoverSeed = trajSeedBase + 500000 + 10000*drop;

        % Select the three distinct SU directions ONCE for this spatial drop.

        % For later T values, keep the same SU identities so the flight-time

        % sweep changes only the available mission duration, not the mission

        % targets themselves. The reachable service point along each retained

        % SU direction is recomputed for the current N below.

        if t_idx == 1

            [q_hover_all,targetSUIdx] = choose_joint_distinct_su_hover_points_3uav( ...
                q_start_all,WS,WP,WW,w_b,w_wifi,H, ...
                f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
                P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
                GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
                Bsc_lic,Bsc_un,area,hoverSeed,trajEvalSamples, ...
                false,useShadow,useShadowA2G,N,STEP_MAX);

            targetSUIdx_fixed = targetSUIdx;

        else

            targetSUIdx = targetSUIdx_fixed;

            q_hover_all = recompute_fixed_su_service_points_3uav( ...
                targetSUIdx,q_start_all,WS,N,STEP_MAX,area);

        end

        % ---------------------------------------------------------------

        % 2) Synchronized SU-oriented non-retracing closed-loop initialization.

        % ---------------------------------------------------------------

        Qhover = zeros(N,3,U);

        for u = 1:U

            Qhover(:,:,u) = build_single_fly_hover_trajectory( ...
                q_start_all(u,:),q_hover_all(u,:),q_end_all(u,:), ...
                N,STEP_MAX,area,u);

        end

        % Jointly shrink only as much as needed if independently generated

        % SU-oriented loops would otherwise approach one another too closely.

        % This operation preserves each UAV's exact start=end condition, loop

        % shape, SU direction, altitude box and speed feasibility. It is used

        % only for initialization; subsequent motion is optimized by SCA.

        Qhover = make_joint_initialization_collision_safe( ...
            Qhover,q_start_all,area,STEP_MAX,UAV_DMIN);

        % ---------------------------------------------------------------

        % 3) Previous-T warm start + proper common-random-number mapping.

        %

        % For a longer flight time, preserve the previous optimized mission and

        % insert the additional time samples by hovering at the PREVIOUS

        % trajectory's best DETERMINISTIC large-scale joint waypoint.  This is

        % used only to construct a good feasible warm start; it does NOT impose

        % any monotonicity requirement on the final Monte-Carlo rate curve.

        %

        % Old slots retain their old RNG identities (common random numbers),

        % while every newly inserted slot receives a fresh unique RNG identity.

        % Therefore no favorable fading realization is copied/duplicated.

        % ---------------------------------------------------------------

        if isempty(Qprev)

            Qinit = Qhover;

            initText = 'SU-oriented closed-loop';

            slotSeedCurrent = 1:N;

        else

            Nold = size(Qprev,1);

            if N < Nold

                error('T_list must be nondecreasing for nested-T warm starts.');

            end

            initSeed = trajSeedBase + 700000 + 10000*drop;

            % Select the insertion point using ONLY deterministic large-scale

            % CSI. Fast fading / Monte-Carlo outcomes are deliberately excluded

            % from this warm-start decision.

            bestJointIdx = find_best_joint_waypoint_large_scale_3uav( ...
                Qprev,WS,WP,WW,w_b,w_wifi,H, ...
                f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
                P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
                GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
                Bsc_lic,Bsc_un,initSeed,trajEvalSamples,useShadow,useShadowA2G);

            Qnested = expand_joint_trajectory_with_hover( ...
                Qprev,N,bestJointIdx);

            % Preserve the RNG identities of inherited old slots, but assign

            % fresh unique identities to newly inserted slots.  This provides

            % paired/common random numbers without cloning a favorable fade.

            slotSeedCurrent = expand_slot_seed_ids( ...
                slotSeedPrev,N,bestJointIdx);

            % Do NOT project or repair Qnested here. The nested construction

            % deliberately keeps every old waypoint exactly and only duplicates

            % one already-feasible synchronized joint waypoint. Hence it

            % preserves the previous mission's speed, regions, endpoints and

            % continuous collision safety without changing its objective.

            for u=1:U

                st = sqrt(sum(diff(Qnested(:,:,u),1,1).^2,2));

                if max(st) > STEP_MAX*(1+1e-8)

                    error('Nested previous-T trajectory violates speed for UAV %d.',u);

                end

                if any(Qnested(:,1,u)<UAV_XMIN(u)-1e-8) || ...
                   any(Qnested(:,1,u)>UAV_XMAX(u)+1e-8)

                    error('Nested previous-T trajectory left the horizontal flight area for UAV %d.',u);

                end

                if norm(Qnested(1,:,u)-q_start_all(u,:))>1e-8 || ...
                   norm(Qnested(end,:,u)-q_end_all(u,:))>1e-8

                    error('Nested previous-T trajectory changed endpoint for UAV %d.',u);

                end

            end

            if min_joint_separation_3uav(Qnested) < UAV_DMIN*(1-1e-8) || ...
               min_continuous_segment_separation_3uav(Qnested) < UAV_DMIN*(1-1e-8) || ...
               min_continuous_segment_horizontal_separation_3uav(Qnested) < UAV_DMIN_XY*(1-1e-8)

                error('Nested previous-T trajectory lost inter-UAV separation.');

            end

            % Use the nested previous-T optimized mission as the warm start.

            % It is already physically feasible and preserves all inherited

            % waypoints. The fresh SU-oriented loop is retained only as the

            % first-T initializer, not as a competing geometry at every T.

            Qinit = Qnested;

            initText = sprintf('previous-T warm start + deterministic hover@%d',bestJointIdx);

        end

        fprintf('\nDrop %d/%d | T=%3d s | N=%3d | init=%s\n', ...
            drop,NDrops,T_fixed,N,initText);

        % ---------------------------------------------------------------

        % 4) Joint analytical Taylor-SCA trajectory optimization.

        %    Resources are fixed inside each local SCA step and refreshed at

        %    the next SCA iteration (alternating-optimization structure).

        % ---------------------------------------------------------------

        [Qprop,convHistory,Qiter] = run_joint_traj_quiet_3uav( ...
            Qinit,q_start_all,q_end_all,WS,WP,WW,w_b,w_wifi,H, ...
            f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
            P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
            GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
            Bsc_lic,Bsc_un,AO_traj_iter,STEP_MAX,alpha_traj,area, ...
            trajSeedBase + 10000*drop + 1000*t_idx, ...
            trajEvalSamples,useRayleigh,useShadow,useShadowA2G);

        convTraceDrop(drop,:,t_idx) = convHistory(:).';
        fprintf('  Final deterministic AO rate = %.6f Mbps\n',convHistory(end));
        fprintf('  Initial deterministic AO rate = %.6f Mbps\n',convHistory(1));

        % TRUE MC convergence: independently evaluate each accepted iterate.
        % Qiter(:,:, :,r) is the actual trajectory after AO/SCA iteration r.
        for r=1:AO_traj_iter
            if r>1 && isequal(Qiter(:,:,:,r),Qiter(:,:,:,r-1))
                rateMC = lastRateMC; % Identical trajectory and random realizations.
            else
            [rateMC,~,~,~,~] = evaluate_final_mc_trajectory_3uav( ...
                Qiter(:,:,:,r),slotSeedCurrent,MC,baseSeed,WS,WP,WW,w_b,w_wifi,H, ...
                f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
                P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
                GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un, ...
                useRayleigh,useShadow,useShadowA2G);
            end
            lastRateMC = rateMC;
            accRateDrop(r) = accRateDrop(r)+rateMC;
            fprintf('  Drop %d/%d, AO iter %d/%d, MC rate = %.6f Mbps\n', ...
                drop,NDrops,r,AO_traj_iter,rateMC);
        end
        % Persist optimized mission and inherited fading identities for next T.
        Qprev = Qprop;
        slotSeedPrev = slotSeedCurrent;
    end

    % Per-iteration MC values are already accumulated above.
    accTxPowerDrop = accTxPowerDrop + avgTxPower;
    accActiveDrop = accActiveDrop + avgActive;
    accModeFracDrop = accModeFracDrop + modeFrac;
    fprintf('Finished flight-time sweep drop %d/%d\n',drop,NDrops);
end

meanRate = accRateDrop/NDrops;
proposedRates = meanRate(:);
end
