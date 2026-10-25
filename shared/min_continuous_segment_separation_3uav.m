function dminNow = min_continuous_segment_separation_3uav(Q)

% Exact closest-approach distance for synchronized, piecewise-linear motion.

% During interval n -> n+1 each UAV is assumed to move linearly using the

% same normalized time tau in [0,1]. For each UAV pair,

%   r(tau) = r0 + tau*dr,

% and the minimizing tau has a closed form and is clipped to [0,1].

[N,~,U] = size(Q);

dminNow = inf;

% N=1 is not used by this simulation, but keep the function well-defined.

if N < 2

    dminNow = min_joint_separation_3uav(Q);

    return;

end

for n=1:N-1

    for u=1:U-1

        for v=u+1:U

            r0 = reshape(Q(n,:,u)   - Q(n,:,v),1,3);

            du = reshape(Q(n+1,:,u) - Q(n,:,u),1,3);

            dv = reshape(Q(n+1,:,v) - Q(n,:,v),1,3);

            dr = du - dv;

            denom = dot(dr,dr);

            if denom > 1e-14

                tauStar = -dot(r0,dr)/denom;

                tauStar = min(max(tauStar,0),1);

            else

                tauStar = 0;

            end

            dNow = norm(r0 + tauStar*dr);

            dminNow = min(dminNow,dNow);

        end

    end

end

end
