function [Q,ok] = make_broad_service_loop(qStart,qServ,N,stepMax,area,u,halfSpan)

% Construct a non-retracing, broad, piecewise-smooth closed route around the

% selected SU-side point.  The route is deliberately expanded along y because

% both x and y use the common flight area.

global UAV_HMIN UAV_HMAX UAV_XMIN UAV_XMAX

qStart=qStart(:).'; qServ=qServ(:).';

% Shift the service-loop centre if needed so both y anchors fit in [0,area].

yc=min(max(qServ(2),halfSpan),area-halfSpan);

yLo=yc-halfSpan; yHi=yc+halfSpan;

% Transit x lies between start and service, keeping both outbound/return paths

% inside the allowed flight area while producing two visibly distinct branches.

xTransit=0.45*qStart(1)+0.55*qServ(1);

xTransit=min(max(xTransit,UAV_XMIN(u)),UAV_XMAX(u));

zMid=0.5*(qStart(3)+qServ(3));

P=[ ...
    qStart; ...
    xTransit, yLo, zMid; ...
    qServ(1), yLo, qServ(3); ...
    qServ; ...
    qServ(1), yHi, qServ(3); ...
    xTransit, yHi, zMid; ...
    qStart];

% Clip only for numerical safety; construction already respects x/y boxes.

for ii=1:size(P,1), P(ii,:)=clip_q3d_role(P(ii,:),area,u); end

% Remove duplicate consecutive anchors that may occur near y boundaries.

d=sqrt(sum(diff(P,1,1).^2,2));

keep=[true; d>1e-9]; P=P(keep,:);

if size(P,1)<3, Q=[]; ok=false; return; end

seg=sqrt(sum(diff(P,1,1).^2,2));

L=[0;cumsum(seg)]; totalL=L(end);

if totalL>(N-1)*stepMax*(1+1e-10)

    Q=[]; ok=false; return;

end

% Arc-length resampling.  'pchip' gives a visually smooth path without the

% large overshoot of a free cubic spline.  If pchip leaves the feasible box,

% fall back to linear interpolation, which is guaranteed to stay in the convex

% region box.

sq=linspace(0,totalL,N).';

Q=zeros(N,3);

for dd=1:3, Q(:,dd)=interp1(L,P(:,dd),sq,'pchip'); end

Q(1,:)=qStart; Q(end,:)=qStart;

inside=all(Q(:,1)>=UAV_XMIN(u)-1e-9 & Q(:,1)<=UAV_XMAX(u)+1e-9 & ...
           Q(:,2)>=-1e-9 & Q(:,2)<=area+1e-9 & ...
           Q(:,3)>=UAV_HMIN-1e-9 & Q(:,3)<=UAV_HMAX+1e-9);

if ~inside

    for dd=1:3, Q(:,dd)=interp1(L,P(:,dd),sq,'linear'); end

    Q(1,:)=qStart; Q(end,:)=qStart;

end

steps=sqrt(sum(diff(Q,1,1).^2,2));

inside=all(Q(:,1)>=UAV_XMIN(u)-1e-9 & Q(:,1)<=UAV_XMAX(u)+1e-9 & ...
           Q(:,2)>=-1e-9 & Q(:,2)<=area+1e-9 & ...
           Q(:,3)>=UAV_HMIN-1e-9 & Q(:,3)<=UAV_HMAX+1e-9);

ok=inside && max(steps)<=stepMax*(1+1e-8) && norm(Q(1,:)-Q(end,:))<=1e-9;

if ~ok, Q=[]; end

end
