function [bestRate,bestMode,bestInfo] = eval_slot_3uav( ...
        Qslot,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,useRayleigh,useShadow,useShadowA2G)

% -------------------------------------------------------------------------

% Instantaneous-CSI mode-control logic:

%   1) Generate the physical channels ONCE for this slot. If useRayleigh=true,

%      each link uses G = beta*|g|^2 with an independent Rayleigh coefficient.

%   2) Different UAVs and the two frequency bands have different channels.

%   3) Reuse these already-generated physical-link channels consistently when

%      evaluating C/C, L/U, C/U and L/C. No new fading is drawn per mode.

%   4) Run the resource allocator for all four modes on that common slot CSI and

%      select the feasible mode with the largest instantaneous sum rate.

%

% Trajectory optimization is intentionally separate: its SCA block uses

% deterministic large-scale Al-Hourani CSI because UAV motion is slow relative

% to fast Rayleigh fading.

% -------------------------------------------------------------------------

% Generate one instantaneous, band-specific channel state for every physical

% link in this slot. This single set is shared by all candidate-mode tests.

[huL,huU,hp2,hw2,hup,huw] = gen_all_channels_3uav( ...
    Qslot,WS,WP,WW,w_b,w_wifi,H, ...
    f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
    useRayleigh,useShadow,useShadowA2G);

% Evaluate all four candidate modes using exactly the same physical-link CSI.

[modeRates,modeInfo] = evaluate_all_modes_from_channels( ...
    huL,huU,hp2,hw2,hup,huw, ...
    P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
    GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un);

% Instantaneous-CSI adaptive mode selection.

[bestRate,bestMode] = max(modeRates);

bestInfo = modeInfo(bestMode);

end
