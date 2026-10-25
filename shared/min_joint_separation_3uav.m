function dminNow = min_joint_separation_3uav(Q)

[N,~,U]=size(Q); dminNow=inf;

for n=1:N

    for u=1:U-1

        for v=u+1:U, dminNow=min(dminNow,norm(Q(n,:,u)-Q(n,:,v))); end

    end

end

end
