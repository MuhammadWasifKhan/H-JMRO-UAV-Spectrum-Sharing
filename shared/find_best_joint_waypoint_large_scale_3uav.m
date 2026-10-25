function bestIdx = find_best_joint_waypoint_large_scale_3uav( ...
        Q,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,seedBase,nAvg,useShadow,useShadowA2G)

% Select the synchronized waypoint of the PREVIOUS optimized trajectory with

% the largest deterministic large-scale post-allocation system rate. The same

% joint waypoint is duplicated for all UAVs when extra T slots are introduced.

Nold = size(Q,1);

bestVal = -inf;

bestIdx = max(1,round(Nold/2));

for n=1:Nold

    Qslot = squeeze(Q(n,:,:)).';

    val = eval_slot_objective_3uav_avg( ...
        Qslot,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,seedBase+n,nAvg, ...
        false,useShadow,useShadowA2G);

    if val > bestVal

        bestVal = val;

        bestIdx = n;

    end

end

end
