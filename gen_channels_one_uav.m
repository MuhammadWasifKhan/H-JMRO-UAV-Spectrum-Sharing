function [huL,huU,hp2,hw2,hup,huw] = gen_channels_one_uav( ...
        q,u,WS,WP,WW,w_b,w_wifi,H, ...
        f_lic,f_un,f_pbs,f_wifi,c,A2G,n_pbs,n_wifi,Shadow, ...
        useRayleigh,useShadow,useShadowA2G)

K=size(WS,1); J=size(WP,1); M=size(WW,1); huL=zeros(K,J); huU=zeros(K,M); hup=zeros(J,1); huw=zeros(M,1);

for k=1:K

    gL=a2g_gain_alhourani(q,WS(k,:),H,f_lic,c,A2G,useShadowA2G,Shadow.a2g_su_z(k,u));

    gU=a2g_gain_alhourani(q,WS(k,:),H,f_un,c,A2G,useShadowA2G,Shadow.a2g_su_z(k,u));

    for j=1:J, if useRayleigh, z=(randn+1j*randn)/sqrt(2); huL(k,j)=gL*abs(z)^2; else, huL(k,j)=gL; end, end

    for m=1:M, if useRayleigh, z=(randn+1j*randn)/sqrt(2); huU(k,m)=gU*abs(z)^2; else, huU(k,m)=gU; end, end

end

for j=1:J

    g=a2g_gain_alhourani(q,WP(j,:),H,f_lic,c,A2G,useShadowA2G,Shadow.a2g_pu_z(j,u)); if useRayleigh, z=(randn+1j*randn)/sqrt(2); hup(j)=g*abs(z)^2; else, hup(j)=g; end

end

for m=1:M

    g=a2g_gain_alhourani(q,WW(m,:),H,f_un,c,A2G,useShadowA2G,Shadow.a2g_wu_z(m,u)); if useRayleigh, z=(randn+1j*randn)/sqrt(2); huw(m)=g*abs(z)^2; else, huw(m)=g; end

end

hp2=zeros(K,J); PL1=20*log10(4*pi*f_pbs/c);

for k=1:K

    d=max(norm(w_b-WS(k,:)),1.0); X=0; if useShadow, X=Shadow.pbs_su_dB(k); end; g=10^(-(PL1+10*n_pbs*log10(d)+X)/10);

    for j=1:J, if useRayleigh, z=(randn+1j*randn)/sqrt(2); hp2(k,j)=g*abs(z)^2; else, hp2(k,j)=g; end, end

end

hw2=zeros(K,M); PL1=20*log10(4*pi*f_wifi/c);

for k=1:K

    d=max(norm(w_wifi-WS(k,:)),1.0); X=0; if useShadow, X=Shadow.wifi_su_dB(k); end; g=10^(-(PL1+10*n_wifi*log10(d)+X)/10);

    for m=1:M, if useRayleigh, z=(randn+1j*randn)/sqrt(2); hw2(k,m)=g*abs(z)^2; else, hw2(k,m)=g; end, end

end

end
