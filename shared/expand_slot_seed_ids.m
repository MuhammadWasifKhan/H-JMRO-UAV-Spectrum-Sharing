function seedNew = expand_slot_seed_ids(seedOld,Nnew,bestIdx)

% Preserve RNG identities of inherited old slots across T (proper CRN pairing).

% Every NEW inserted slot receives a fresh unique RNG identity; no previous

% fading realization is duplicated merely because its geometry is duplicated.

seedOld = seedOld(:).';

Nold = numel(seedOld);

if Nnew < Nold

    error('Nnew must be >= numel(seedOld).');

end

if Nold==0

    seedNew = 1:Nnew;

    return;

end

if Nnew==Nold

    seedNew = seedOld;

    return;

end

extra = Nnew-Nold;

bestIdx = min(max(round(bestIdx),1),Nold);

nextId = max(seedOld) + 1;

freshIds = nextId:(nextId+extra-1);

seedNew = [seedOld(1:bestIdx), freshIds, seedOld(bestIdx+1:end)];

end
