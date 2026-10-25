function dminNow = min_full_path_horizontal_separation_3uav(Q)

% Minimum geometric clearance between the COMPLETE XY path traces of every

% UAV pair. Unlike synchronized collision checking, this compares every

% segment of one UAV with every segment of the other UAV, irrespective of

% slot/time index. Hence two plotted trajectories cannot cross at different

% times and still be declared safe.

[N,~,U] = size(Q);

dminNow = inf;

if N < 2

    for u=1:U-1

        for v=u+1:U

            dminNow = min(dminNow,norm(Q(1,1:2,u)-Q(1,1:2,v)));

        end

    end

    return;

end

for u=1:U-1

    for v=u+1:U

        for n=1:N-1

            p1 = reshape(Q(n,1:2,u),1,2);

            p2 = reshape(Q(n+1,1:2,u),1,2);

            for m=1:N-1

                q1 = reshape(Q(m,1:2,v),1,2);

                q2 = reshape(Q(m+1,1:2,v),1,2);

                dminNow = min(dminNow,segment_segment_distance_2d(p1,p2,q1,q2));

                if dminNow <= 1e-12, return; end

            end

        end

    end

end

end
