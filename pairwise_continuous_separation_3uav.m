function Dpair = pairwise_continuous_separation_3uav(Q)

% Exact minimum synchronized segment distance for each UAV pair.

[N,~,U] = size(Q);

Dpair = inf(U,U);

for u=1:U

    Dpair(u,u)=0;

end

if N < 2

    for u=1:U-1

        for v=u+1:U

            d=norm(Q(1,:,u)-Q(1,:,v));

            Dpair(u,v)=d; Dpair(v,u)=d;

        end

    end

    return;

end

for u=1:U-1

    for v=u+1:U

        dm=inf;

        for n=1:N-1

            r0 = reshape(Q(n,:,u)-Q(n,:,v),1,3);

            du = reshape(Q(n+1,:,u)-Q(n,:,u),1,3);

            dv = reshape(Q(n+1,:,v)-Q(n,:,v),1,3);

            dr = du-dv;

            den = dot(dr,dr);

            if den>1e-14

                tau = min(max(-dot(r0,dr)/den,0),1);

            else

                tau = 0;

            end

            dm=min(dm,norm(r0+tau*dr));

        end

        Dpair(u,v)=dm; Dpair(v,u)=dm;

    end

end

end
