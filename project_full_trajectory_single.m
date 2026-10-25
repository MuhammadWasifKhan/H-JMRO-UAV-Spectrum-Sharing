function Q = project_full_trajectory_single( ...
        Q,area,stepMax,q_start,q_end,projIter,u)

N = size(Q,1);

for nn=1:size(Q,1)

    Q(nn,:) = clip_q3d_role(Q(nn,:),area,u);

end

Q(1,:) = q_start;

Q(N,:) = q_end;

for it=1:projIter

    Q(1,:) = q_start;

    for n=2:N

        d = Q(n,:)-Q(n-1,:);

        nd = norm(d);

        if nd>stepMax

            Q(n,:) = Q(n-1,:) + stepMax*d/max(nd,1e-12);

        end

        Q(n,:) = clip_q3d_role(Q(n,:),area,u);

    end

    Q(N,:) = q_end;

    for n=N-1:-1:1

        d = Q(n,:)-Q(n+1,:);

        nd = norm(d);

        if nd>stepMax

            Q(n,:) = Q(n+1,:) + stepMax*d/max(nd,1e-12);

        end

        Q(n,:) = clip_q3d_role(Q(n,:),area,u);

    end

    Q(1,:) = q_start;

    Q(N,:) = q_end;

end

end

% ========================== TRUE-3D clipping helper ======================
