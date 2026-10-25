function q = clip_q3d(q,area)

global UAV_HMIN UAV_HMAX

q = double(q(:).');

q(1) = min(max(q(1),0),area);

q(2) = min(max(q(2),0),area);

q(3) = min(max(q(3),UAV_HMIN),UAV_HMAX);

end
