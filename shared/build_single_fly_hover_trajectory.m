function Q = build_single_fly_hover_trajectory( ...
        q_start,q_hover,q_end,N,stepMax,area,u)

% =========================================================================

% BROAD SU-SERVICE CLOSED LOOP (ICC FINAL)

% =========================================================================

% Earlier ellipse construction could still produce a compact red/yellow loop

% when the selected service point was close to the fixed start point.  Here the

% closed route is built explicitly around the selected SU-side service point:

%

% start -> lower transit -> lower service anchor -> exact service point

%       -> upper service anchor -> upper transit -> start.

%

% The two service anchors are separated mainly along y, where every UAV has

% the full [0,1000] m freedom.  Hence the loop remains visibly expanded even

% when the SU and start have similar x coordinates.  The x coordinate always

% stays inside the role-specific segment.  The route is then arc-length

% resampled, so the inter-waypoint displacement is directly controlled by the

% flight-time/speed budget.  If necessary, only the loop half-width is reduced;

% the selected SU-side service point is retained whenever physically feasible.

% =========================================================================

global UAV_HMIN UAV_HMAX UAV_XMIN UAV_XMAX

q_start = clip_q3d_role(q_start(:).',area,u);

q_hover = clip_q3d_role(q_hover(:).',area,u);

q_end   = clip_q3d_role(q_end(:).',area,u);

if norm(q_end-q_start)>1e-8

    error('Broad service-loop initializer requires q_end=q_start for UAV %d.',u);

end

if N<3

    Q=repmat(q_start,N,1); return;

end

% Exact in-region service point nearest the selected SU.

q_hover(1)=min(max(q_hover(1),UAV_XMIN(u)),UAV_XMAX(u));

q_hover(2)=min(max(q_hover(2),0),area);

q_hover(3)=min(max(q_hover(3),UAV_HMIN),UAV_HMAX);

% Try a broad loop first.  120 m half-span means ~240 m visible y spread.

% Reduce only if the available flight distance cannot support it.

halfSpanList = [120 110 100 90 80 70 60 50 40 30 20]; % MUST match reachability test

bestQ=[]; usedHalfSpan=0; reachedTarget=false;

for hs=halfSpanList

    [Qtry,ok] = make_broad_service_loop(q_start,q_hover,N,stepMax,area,u,hs);

    if ok

        bestQ=Qtry; usedHalfSpan=hs; reachedTarget=true; break;

    end

end

% Reachability-aware target selection now tests the SAME broad-loop family

% used here (minimum half-span 20 m). Therefore any selected target/service

% point is guaranteed to be constructible by this initializer; no zero-width

% out-and-back path can pass the selector.

if isempty(bestQ)

    error(['Internal reachability inconsistency for UAV %d. The selector and ' ...
           'broad-loop initializer should now use identical span tests.'],u);

end

Q=bestQ;

Q(1,:)=q_start; Q(end,:)=q_start;

steps=sqrt(sum(diff(Q,1,1).^2,2));

if max(steps)>stepMax*(1+1e-8), error('Broad service loop violates speed for UAV %d.',u); end

if any(Q(:,1)<UAV_XMIN(u)-1e-9) || any(Q(:,1)>UAV_XMAX(u)+1e-9)

    error('Broad service loop left assigned region for UAV %d.',u);

end

targetGap=min(sqrt(sum((Q(:,1:2)-q_hover(1:2)).^2,2)));

xySpan=[max(Q(:,1))-min(Q(:,1)),max(Q(:,2))-min(Q(:,2))];

fprintf(['    UAV %d broad service loop: target-gap=%.2f m | y-half-span=%.1f m | ' ...
         'XY span=[%.1f, %.1f] m | full-target=%d\n'], ...
         u,targetGap,usedHalfSpan,xySpan(1),xySpan(2),reachedTarget);

end
