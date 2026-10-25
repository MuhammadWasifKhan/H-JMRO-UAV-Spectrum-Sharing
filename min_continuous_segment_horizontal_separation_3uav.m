function dminNow = min_continuous_segment_horizontal_separation_3uav(Q)

% Exact closest HORIZONTAL approach for synchronized piecewise-linear motion.

% This prevents the optimizer from satisfying collision avoidance only by

% putting two UAVs at different altitudes while their x-y paths overlap.

[N,~,U] = size(Q);

dminNow = inf;

if N < 2

    for u=1:U-1

        for v=u+1:U

            dminNow=min(dminNow,norm(Q(1,1:2,u)-Q(1,1:2,v)));

        end

    end

    return;

end

for n=1:N-1

    for u=1:U-1

        for v=u+1:U

            r0 = reshape(Q(n,1:2,u)-Q(n,1:2,v),1,2);

            du = reshape(Q(n+1,1:2,u)-Q(n,1:2,u),1,2);

            dv = reshape(Q(n+1,1:2,v)-Q(n,1:2,v),1,2);

            dr = du-dv;

            den = dot(dr,dr);

            if den > 1e-14

                tau = min(max(-dot(r0,dr)/den,0),1);

            else

                tau = 0;

            end

            dminNow=min(dminNow,norm(r0+tau*dr));

        end

    end

end

end
