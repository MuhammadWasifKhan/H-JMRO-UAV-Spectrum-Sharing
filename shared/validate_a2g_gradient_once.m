function validate_a2g_gradient_once(P,c,f,useShadowA2G)

q=[350,620,110]; w=[520,470]; zc=0.37; epsFD=1e-3;

[g,ga]=a2g_gain_alhourani_grad_3d(q,w,f,c,P,useShadowA2G,zc); gf=zeros(1,3);

for d=1:3

    qp=q; qm=q; qp(d)=qp(d)+epsFD; qm(d)=qm(d)-epsFD;

    gp=a2g_gain_alhourani(qp,w,100,f,c,P,useShadowA2G,zc); gm=a2g_gain_alhourani(qm,w,100,f,c,P,useShadowA2G,zc);

    gf(d)=(gp-gm)/(2*epsFD);

end

rel=norm(ga-gf)/max(norm(gf),1e-20); fprintf('Analytical A2G gradient check: relative error = %.3e\n',rel);

if rel>5e-3, warning('A2G analytical gradient validation error is larger than expected.'); end

if ~isfinite(g), error('Gradient validation produced non-finite channel gain.'); end

end
