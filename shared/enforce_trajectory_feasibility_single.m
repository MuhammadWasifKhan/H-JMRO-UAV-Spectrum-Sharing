function Q = enforce_trajectory_feasibility_single( ...
        Q,area,stepMax,q_start,q_end,u)

% Guaranteed feasibility repair for a whole trajectory.

%

% 1) Clip to the 3-D flight box and fix endpoints.

% 2) Try alternating forward/backward projections.

% 3) If a small residual remains, blend the candidate trajectory with the

%    straight feasible start-to-end path. Because the speed-ball and box

%    constraints are convex, the blend is feasible for a sufficiently small

%    coefficient and preserves as much of the optimized shape as possible.

N = size(Q,1);

q_start = clip_q3d_role(q_start,area,u);

q_end   = clip_q3d_role(q_end,area,u);

% Basic mission feasibility check.

if norm(q_end-q_start) > stepMax*(N-1) + 1e-9

    error(['Mission endpoints are infeasible: start-to-end distance %.6f m ' ...
           'exceeds available travel distance %.6f m.'], ...
           norm(q_end-q_start),stepMax*(N-1));

end

for n=1:N

    Q(n,:) = clip_q3d_role(Q(n,:),area,u);

end

Q(1,:) = q_start;

Q(end,:) = q_end;

% Fast alternating repair.

Q = project_full_trajectory_single(Q,area,stepMax,q_start,q_end,200,u);

steps = sqrt(sum(diff(Q,1,1).^2,2));

if isempty(steps) || max(steps) <= stepMax + 1e-8

    Q(1,:) = q_start;

    Q(end,:) = q_end;

    return;

end

% Straight start-to-end reference; always feasible when the mission itself is.

t = linspace(0,1,N).';

Qlin = (1-t).*repmat(q_start,N,1) + t.*repmat(q_end,N,1);

% Bisection for the largest shape-retention coefficient beta.

lo = 0;

hi = 1;

for it=1:70

    beta = 0.5*(lo+hi);

    Qt = Qlin + beta*(Q-Qlin);

    Qt(1,:) = q_start;

    Qt(end,:) = q_end;

    st = sqrt(sum(diff(Qt,1,1).^2,2));

    if max(st) <= stepMax + 1e-10

        lo = beta;

    else

        hi = beta;

    end

end

Q = Qlin + lo*(Q-Qlin);

Q(1,:) = q_start;

Q(end,:) = q_end;

for n=1:N

    Q(n,:) = clip_q3d_role(Q(n,:),area,u);

end

% Final conservative pass and check.

Q = project_full_trajectory_single(Q,area,stepMax,q_start,q_end,50,u);

steps = sqrt(sum(diff(Q,1,1).^2,2));

if ~isempty(steps) && max(steps) > stepMax + 1e-7

    % Last-resort feasible path. This should never be needed, but it is much

    % safer than terminating a long Monte-Carlo run because of solver residue.

    Q = Qlin;

end

Q(1,:) = q_start;

Q(end,:) = q_end;

end
