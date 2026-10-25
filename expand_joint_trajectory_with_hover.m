function Qnew = expand_joint_trajectory_with_hover(Qold,Nnew,bestIdx)

% Nested-T construction: keep the complete old synchronized mission exactly

% and insert the additional samples by hovering at one existing joint waypoint.

[Nold,D,U] = size(Qold);

if Nnew < Nold

    error('Nnew must be >= Nold.');

end

if Nnew == Nold

    Qnew = Qold;

    return;

end

extra = Nnew-Nold;

bestIdx = min(max(round(bestIdx),1),Nold);

Qnew = zeros(Nnew,D,U);

Qnew(1:bestIdx,:,:) = Qold(1:bestIdx,:,:);

for e=1:extra

    Qnew(bestIdx+e,:,:) = Qold(bestIdx,:,:);

end

if bestIdx < Nold

    Qnew(bestIdx+extra+1:end,:,:) = Qold(bestIdx+1:end,:,:);

end

end
