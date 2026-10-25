function q = feasible_segment_from_q0_single( ...
        q0,qtarget,q_prev,q_next_fixed,stepMax,area)

% Return the furthest point on the segment q0 -> qtarget that satisfies the

% two neighbouring speed constraints and the 3-D box. q0 is assumed feasible.

q0 = clip_q3d(q0,area);

qtarget = clip_q3d(qtarget,area);

q_prev = double(q_prev(:).');

q_next_fixed = double(q_next_fixed(:).');

% If the solver point is already feasible, keep it.

if norm(qtarget-q_prev) <= stepMax + 1e-9 && ...
   norm(qtarget-q_next_fixed) <= stepMax + 1e-9

    q = qtarget;

    return;

end

% q0 should be feasible. If round-off made it marginally infeasible, repair it

% first with the conservative intersection projection.

if norm(q0-q_prev) > stepMax + 1e-9 || ...
   norm(q0-q_next_fixed) > stepMax + 1e-9

    q0 = project_between_neighbors_single(q0,q_prev,q_next_fixed,stepMax,area);

end

lo = 0;

hi = 1;

for it = 1:60

    mid = 0.5*(lo+hi);

    qcand = q0 + mid*(qtarget-q0);

    qcand = clip_q3d(qcand,area);

    if norm(qcand-q_prev) <= stepMax + 1e-10 && ...
       norm(qcand-q_next_fixed) <= stepMax + 1e-10

        lo = mid;

    else

        hi = mid;

    end

end

q = q0 + lo*(qtarget-q0);

q = clip_q3d(q,area);

end
