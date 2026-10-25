function q = clip_q3d_role(q,area,u)

% Clip a UAV waypoint to its allowed 3-D operating region.

% In the current ICC configuration UAV_XMIN=0 and UAV_XMAX=area for all UAVs.

global UAV_HMIN UAV_HMAX UAV_XMIN UAV_XMAX

q = double(q(:).');

q(1) = min(max(q(1),UAV_XMIN(u)),UAV_XMAX(u));

q(2) = min(max(q(2),0),area);

q(3) = min(max(q(3),UAV_HMIN),UAV_HMAX);

end

% =================== One-to-one Hungarian assignment =====================
