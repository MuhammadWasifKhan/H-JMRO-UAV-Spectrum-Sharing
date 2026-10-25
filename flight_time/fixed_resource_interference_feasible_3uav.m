function ok = fixed_resource_interference_feasible_3uav( ...
        Q,fixedInfo,WP,WW,f_lic,f_un,c,A2G,Shadow, ...
        GammaLic,GammaUn,useShadowA2G)

% Check the ORIGINAL nonlinear individual cochannel interference limits while keeping

% the SCA iteration's resource variables fixed. This is the feasibility guard

% used after solving the affine Taylor subproblem.

[N,~,~] = size(Q);

ok = true;

for n=1:N

    info = fixedInfo{n};

    IL = 0;

    u = info.txL;

    for ell=1:size(info.idxLic,1)

        ss = info.idxLic(ell,2);

        p  = info.pLic(ell);

        [g,~] = a2g_gain_alhourani_grad_3d( ...
            Q(n,:,u),WP(ss,:),f_lic,c,A2G,useShadowA2G,Shadow.a2g_pu_z(ss,u));

        IL = max(IL,p*g);

    end

    IU = 0;

    u = info.txU;

    for ell=1:size(info.idxUn,1)

        ss = info.idxUn(ell,2);

        p  = info.pUn(ell);

        [g,~] = a2g_gain_alhourani_grad_3d( ...
            Q(n,:,u),WW(ss,:),f_un,c,A2G,useShadowA2G,Shadow.a2g_wu_z(ss,u));

        IU = max(IU,p*g);

    end

    if IL > GammaLic*(1+1e-7) || IU > GammaUn*(1+1e-7)

        ok = false;

        return;

    end

end

end
