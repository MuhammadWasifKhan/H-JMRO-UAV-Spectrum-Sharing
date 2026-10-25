function Qsafe = make_joint_initialization_collision_safe(Q,qStart,area,stepMax,dminReq)

% Keep the independently generated SU-oriented closed loops whenever they are

% already safe. If they are not, shrink all intermediate displacements toward

% each UAV's own start/end point by the largest tested common scale that gives

% exact continuous-time d_min safety. Scaling around each distinct start point

% preserves loop closure and direction toward the selected SU service region.

global UAV_DMIN_XY

[N,~,U] = size(Q);

Qsafe = Q;

if min_continuous_segment_separation_3uav(Qsafe) >= dminReq*(1+1e-6) && ...
   min_continuous_segment_horizontal_separation_3uav(Qsafe) >= UAV_DMIN_XY*(1+1e-6)

    return;

end

% Search from almost-full loops downward, retaining the largest safe scale.

scales = [linspace(0.98,0.10,89), 0.08,0.06,0.04,0.02,0.01];

found = false;

for ss = scales

    Qt = zeros(size(Q));

    for u=1:U

        base = reshape(qStart(u,:),1,3);

        Qt(:,:,u) = repmat(base,N,1) + ss*(Q(:,:,u)-repmat(base,N,1));

        Qt(1,:,u)   = base;

        Qt(end,:,u) = base;

        Qt(:,:,u) = enforce_trajectory_feasibility_single( ...
            Qt(:,:,u),area,stepMax,base,base,u);

        Qt(1,:,u)   = base;

        Qt(end,:,u) = base;

    end

    if min_continuous_segment_separation_3uav(Qt) >= dminReq*(1+1e-6) && ...
       min_continuous_segment_horizontal_separation_3uav(Qt) >= UAV_DMIN_XY*(1+1e-6)

        Qsafe = Qt;

        found = true;

        break;

    end

end

if ~found

    error(['Unable to construct a continuously non-colliding SU-oriented ' ...
           'closed-loop initialization. Check the three start positions.']);

end

end
