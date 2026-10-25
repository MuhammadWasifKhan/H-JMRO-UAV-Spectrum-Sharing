function [g_lin,grad_g] = a2g_gain_alhourani_grad_3d( ...
        q,w,f,c,P,useShadowA2G,z_coherent)

% q=[x,y,z], ground terminal w=[xg,yg].

q = double(q(:).');

w = double(w(:).');

dx = q(1)-w(1);

dy = q(2)-w(2);

z  = max(q(3),1e-6);

r2 = dx^2 + dy^2;

r  = sqrt(max(r2,1e-18));

d2 = r2 + z^2;

d  = sqrt(max(d2,1e-18));

theta = atan2d(z,r);

P_LoS  = 1/(1 + P.a*exp(-P.b*(theta-P.a)));

P_NLoS = 1-P_LoS;

PL_FS_dB = 20*log10(4*pi*f*d/c);

eta_mean_dB = P_LoS*P.etaLoS_dB + P_NLoS*P.etaNLoS_dB;

sigma_mix = sqrt(max( ...
    P_LoS*(P.sigmaLoS_dB^2) + ...
    P_NLoS*(P.sigmaNLoS_dB^2),1e-18));

Xsigma_dB = 0;

if useShadowA2G

    Xsigma_dB = sigma_mix*z_coherent;

end

PL_dB = PL_FS_dB + eta_mean_dB + Xsigma_dB;

g_lin = 10^(-PL_dB/10);

% ---- analytical gradient of PL_FS ----

grad_PLFS = (20/log(10)) * [dx,dy,z] / max(d2,1e-18);

% ---- analytical gradient of elevation angle in DEGREES ----

% theta = atan2(z,r)

if r > 1e-9

    grad_theta_rad = [ ...
        -z*dx/(r*max(d2,1e-18)), ...
        -z*dy/(r*max(d2,1e-18)), ...
         r/max(d2,1e-18)];

else

    grad_theta_rad = [0,0,0];

end

grad_theta_deg = (180/pi)*grad_theta_rad;

% ---- derivative of LoS probability w.r.t. theta (theta in degrees) ----

dP_dtheta = P.b*P_LoS*(1-P_LoS);

% Mean excess-loss derivative.

deta_dtheta = (P.etaLoS_dB-P.etaNLoS_dB)*dP_dtheta;

% Shadow standard-deviation derivative.

dsigma_dtheta = 0;

if useShadowA2G

    dsigma_dtheta = ...
        0.5*(P.sigmaLoS_dB^2-P.sigmaNLoS_dB^2) * ...
        dP_dtheta / max(sigma_mix,1e-18);

end

dPL_dtheta = deta_dtheta + z_coherent*dsigma_dtheta;

grad_PL = grad_PLFS + dPL_dtheta*grad_theta_deg;

% g = 10^(-PL/10)

grad_g = -(log(10)/10)*g_lin*grad_PL;

if ~all(isfinite(grad_g))

    grad_g = zeros(1,3);

end

end
