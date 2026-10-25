function Q = repair_joint_separation_3uav(Q,area,stepMax,qStart,qEnd,dminReq)

[N,~,U]=size(Q);

for outer=1:40

    changed=false;

    for n=2:N-1

        for u=1:U-1

            for v=u+1:U

                dvec=Q(n,:,u)-Q(n,:,v); d=norm(dvec);

                if d<dminReq

                    changed=true; if d<1e-9, dir=[1,0,0]; else, dir=dvec/d; end; push=0.5*(dminReq-d+1e-3);

                    Q(n,:,u)=clip_q3d_role(Q(n,:,u)+push*dir,area,u); Q(n,:,v)=clip_q3d_role(Q(n,:,v)-push*dir,area,v);

                end

            end

        end

    end

    for u=1:U, Q(:,:,u)=enforce_trajectory_feasibility_single(Q(:,:,u),area,stepMax,qStart(u,:),qEnd(u,:),u); end

    if ~changed || min_joint_separation_3uav(Q)>=dminReq*(1-1e-6), break; end

end

if min_joint_separation_3uav(Q)<dminReq*(1-1e-4), warning('Initialization remains slightly below waypoint d_min; SCA may be infeasible.'); end

end
