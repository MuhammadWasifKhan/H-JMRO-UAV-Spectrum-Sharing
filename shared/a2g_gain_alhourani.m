function g_lin = a2g_gain_alhourani( ...
        q,w,H,f,c,P,useShadowA2G,z_coherent)

% TRUE-3D A2G geometry. Ground node w=[x,y], UAV q=[x,y,z].

q = q(:).';

w = w(:).';

qxy = q(1:2);

zUAV = q(3);

r = norm(qxy-w);

r = max(r,1e-6);

d = sqrt(zUAV^2+r^2);

theta = atan2d(zUAV,r);

P_LoS = 1/(1+P.a*exp(-P.b*(theta-P.a)));

P_NLoS = 1-P_LoS;

PL_FS_dB = 20*log10(4*pi*f*d/c);

eta_mean_dB = ...
    P_LoS*P.etaLoS_dB + P_NLoS*P.etaNLoS_dB;

Xsigma_dB = 0;

if useShadowA2G

    sigma_mix = sqrt( ...
        P_LoS*(P.sigmaLoS_dB^2) + ...
        P_NLoS*(P.sigmaNLoS_dB^2));

    Xsigma_dB = sigma_mix*z_coherent;

end

PL_dB = PL_FS_dB + eta_mean_dB + Xsigma_dB;

g_lin = 10^(-PL_dB/10);

end
