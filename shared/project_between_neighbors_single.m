function q = project_between_neighbors_single(q,q_prev,q_next,stepMax,area)

% Numerical projection onto the intersection of the two mobility balls and the

% 3D box. For a feasible current trajectory this is normally a no-op.

q = clip_q3d(double(q(:).'),area);

q_prev = double(q_prev(:).');

q_next = double(q_next(:).');

for it=1:30

    d = q-q_prev;

    nd = norm(d);

    if nd>stepMax

        q = q_prev + stepMax*d/max(nd,1e-12);

    end

    q = clip_q3d(q,area);

    d = q-q_next;

    nd = norm(d);

    if nd>stepMax

        q = q_next + stepMax*d/max(nd,1e-12);

    end

    q = clip_q3d(q,area);

end

end
