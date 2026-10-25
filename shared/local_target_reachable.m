function [ok,zBest,spanBest] = local_target_reachable( ...
        q0,suXY,alpha,altCand,Ntraj,stepMax,area,u,hsTry)

% Test a target at q0_xy + alpha*(SU-q0_xy) using the SAME broad-loop

% initializer used by the main trajectory construction.

ok = false;

zBest = NaN;

spanBest = 0;

qXY = q0(1:2) + alpha*(suXY-q0(1:2));

for aa=1:numel(altCand)

    qc=[qXY,altCand(aa)];

    for hs=hsTry

        % IMPORTANT: use exactly the same broad-loop family as the actual

        % initializer. Do not accept a zero-width out-and-back trajectory,

        % because it is not a visible closed loop and the downstream

        % initializer/non-retracing constraints cannot use it.

        [~,okLoop]=make_broad_service_loop(q0,qc,Ntraj,stepMax,area,u,hs);

        if okLoop

            ok=true; zBest=qc(3); spanBest=hs; return;

        end

    end

end

end
