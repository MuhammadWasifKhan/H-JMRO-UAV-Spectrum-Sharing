function ok = fixed_resource_interference_feasible_3uav( ...
        Q,fixedInfo,WP,WW,f_lic,f_un,c,A2G,Shadow, ...
        GammaLic,GammaUn,useShadowA2G)

% Check ORIGINAL nonlinear individual cochannel interference constraints while keeping

% the SCA iteration's resource variables fixed. This is the feasibility guard

% used after solving the affine Taylor subproblem.

[N,~,~] = size(Q);

ok = true;

for n=1:N

    info = fixedInfo{n};

    u = info.txL;
    for ell=1:size(info.idxLic,1)
        ss=info.idxLic(ell,2);
        [g,~]=a2g_gain_alhourani_grad_3d( ...
            Q(n,:,u),WP(ss,:),f_lic,c,A2G,useShadowA2G,Shadow.a2g_pu_z(ss,u));
        if info.pLic(ell)*g>GammaLic*(1+1e-7), ok=false; return; end
    end
    u = info.txU;
    for ell=1:size(info.idxUn,1)
        ss=info.idxUn(ell,2);
        [g,~]=a2g_gain_alhourani_grad_3d( ...
            Q(n,:,u),WW(ss,:),f_un,c,A2G,useShadowA2G,Shadow.a2g_wu_z(ss,u));
        if info.pUn(ell)*g>GammaUn*(1+1e-7), ok=false; return; end
    end
end

end
