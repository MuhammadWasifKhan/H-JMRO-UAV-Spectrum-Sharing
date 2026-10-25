function [D,bestInfo] = build_reference_sca_data_3uav( ...
        Qslot,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
        GammaLic,GammaUn,Pmax_lic,Pmax_un, ...
        Bsc_lic,Bsc_un,useShadow,useShadowA2G)

% Build one-slot data for the EXACT first-order Taylor Al-Hourani subproblem.

% Mode/subchannel assignment/powers are optimized first and then held fixed.

[huL,huU,hp2,hw2,hup,huw] = gen_all_channels_3uav( ...
    Qslot,WS,WP,WW,w_b,w_wifi,H, ...
    f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
    false,useShadow,useShadowA2G);

[modeRates,modeInfo] = evaluate_all_modes_from_channels( ...
    huL,huU,hp2,hw2,hup,huw, ...
    P_PBS_j,P_WIFI_m,sigma2_lic,sigma2_un, ...
    GammaLic,GammaUn,Pmax_lic,Pmax_un,Bsc_lic,Bsc_un);

[~,bestMode] = max(modeRates);

bestInfo = modeInfo(bestMode);

% -------------------------------------------------------------------------

% STRICT-FEASIBILITY HEADROOM FOR TRAJECTORY SCA ONLY

% -------------------------------------------------------------------------

% The KKT resource solution may be active at an individual

% interference boundary. That is optimal for the fixed-position resource

% problem, but it leaves no room for a subsequent trajectory displacement:

% any arbitrarily small increase of the protection-link gain can then make

% the fixed-resource SCA trial infeasible.

%

% Therefore, ONLY while building the fixed-resource reference used by the

% trajectory-SCA subproblem, back off the selected KKT powers just enough to

% leave a small strict-interior interference margin. The normal slot/resource

% optimizer and the final Monte-Carlo evaluation still use the original KKT

% solution and the original Gamma thresholds; nothing is changed there.

trajHeadroom = 0.98;   % keep <=98% of Gamma during trajectory-SCA reference

% Apply 2% strict headroom independently to each selected protected receiver.
for rr=1:size(bestInfo.idxLic,1)
    ss=bestInfo.idxLic(rr,2);
    bound=trajHeadroom*GammaLic/max(hup(ss,bestInfo.txL),realmin);
    bestInfo.pLic(rr)=min(bestInfo.pLic(rr),bound);
end
for rr=1:size(bestInfo.idxUn,1)
    ss=bestInfo.idxUn(rr,2);
    bound=trajHeadroom*GammaUn/max(huw(ss,bestInfo.txU),realmin);
    bestInfo.pUn(rr)=min(bestInfo.pUn(rr),bound);
end

% Refresh power-accounting fields when present. These fields are diagnostic

% only inside the trajectory-SCA reference; mode selection itself is unchanged.

if isfield(bestInfo,'txPowerLic'), bestInfo.txPowerLic = sum(bestInfo.pLic); end

if isfield(bestInfo,'txPowerUn'),  bestInfo.txPowerUn  = sum(bestInfo.pUn);  end

if isfield(bestInfo,'txPowerTotal'), bestInfo.txPowerTotal = sum(bestInfo.pLic)+sum(bestInfo.pUn); end

emptyRate = struct('u',{},'k',{},'ss',{},'p',{},'rate0',{},'gradRate',{});

emptyProt = struct('u',{},'ss',{},'p',{},'I0',{},'gradI',{});

D.lic = emptyRate;

D.un = emptyRate;

D.protLic = emptyProt;

D.protUn = emptyProt;

% ----------------------------- Licensed desired/protection links ----------

u = bestInfo.txL;

for rr=1:size(bestInfo.idxLic,1)

    k  = bestInfo.idxLic(rr,1);

    ss = bestInfo.idxLic(rr,2);

    p  = bestInfo.pLic(rr);

    q0 = Qslot(u,:);

    % Exact Al-Hourani desired-link gain and analytical 3-D gradient.

    [g0,gradG] = a2g_gain_alhourani_grad_3d( ...
        q0,WS(k,:),f_lic,c,A2G,useShadowA2G,Shadow.a2g_su_z(k,u));

    denom = max(hp2(k,ss)*P_PBS_j(ss)+sigma2_lic,1e-30);

    a = max(p/denom,0);

    rate0 = Bsc_lic*log2(1 + a*g0);

    gradRate = (Bsc_lic/log(2)) * (a/max(1+a*g0,1e-30)) * gradG;

    D.lic(end+1) = struct('u',u,'k',k,'ss',ss,'p',p, ...
        'rate0',rate0,'gradRate',gradRate); %#ok<AGROW>

    % Exact Al-Hourani protection-link gain and analytical gradient.

    [gP0,gradGP] = a2g_gain_alhourani_grad_3d( ...
        q0,WP(ss,:),f_lic,c,A2G,useShadowA2G,Shadow.a2g_pu_z(ss,u));

    D.protLic(end+1) = struct('u',u,'ss',ss,'p',p, ...
        'I0',p*gP0,'gradI',p*gradGP); %#ok<AGROW>

end

% ----------------------------- Unlicensed desired/protection links --------

u = bestInfo.txU;

for rr=1:size(bestInfo.idxUn,1)

    k  = bestInfo.idxUn(rr,1);

    ss = bestInfo.idxUn(rr,2);

    p  = bestInfo.pUn(rr);

    q0 = Qslot(u,:);

    [g0,gradG] = a2g_gain_alhourani_grad_3d( ...
        q0,WS(k,:),f_un,c,A2G,useShadowA2G,Shadow.a2g_su_z(k,u));

    denom = max(hw2(k,ss)*P_WIFI_m(ss)+sigma2_un,1e-30);

    a = max(p/denom,0);

    rate0 = Bsc_un*log2(1 + a*g0);

    gradRate = (Bsc_un/log(2)) * (a/max(1+a*g0,1e-30)) * gradG;

    D.un(end+1) = struct('u',u,'k',k,'ss',ss,'p',p, ...
        'rate0',rate0,'gradRate',gradRate); %#ok<AGROW>

    [gW0,gradGW] = a2g_gain_alhourani_grad_3d( ...
        q0,WW(ss,:),f_un,c,A2G,useShadowA2G,Shadow.a2g_wu_z(ss,u));

    D.protUn(end+1) = struct('u',u,'ss',ss,'p',p, ...
        'I0',p*gW0,'gradI',p*gradGW); %#ok<AGROW>

end

end
