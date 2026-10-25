function q = project_single_step_single( ...
        q,q_prev,q_end,n,N,stepMax,area)

q = double(q(:).');

q_prev = double(q_prev(:).');

q_end = double(q_end(:).');

q = clip_q3d(q,area);

d = q-q_prev;

nd = norm(d);

if nd>stepMax

    q = q_prev + stepMax*d/max(nd,1e-12);

end

remSteps = N-n;

rReach = remSteps*stepMax;

v = q-q_end;

nv = norm(v);

if nv>rReach

    q = q_end + rReach*v/max(nv,1e-12);

end

q = clip_q3d(q,area);

d = q-q_prev;

nd = norm(d);

if nd>stepMax

    q = q_prev + stepMax*d/max(nd,1e-12);

end

q = clip_q3d(q,area);

end
