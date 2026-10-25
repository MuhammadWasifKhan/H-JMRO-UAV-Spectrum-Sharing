function pairs = hungarian_one_to_one_max(Umat)

% Maximum-weight one-to-one assignment using the Hungarian algorithm.

%

% INPUT

%   Umat(k,s): utility of assigning SU k to subchannel s.

%

% OUTPUT

%   pairs(:,1): assigned SU indices

%   pairs(:,2): assigned subchannel indices

%

% The routine maximizes the SUM utility subject to

%   - each SU being assigned to at most one subchannel;

%   - each subchannel being assigned to at most one SU.

%

% For the ICC setting K=6 and S=4, all four subchannels are assigned to four

% distinct SUs. The implementation is self-contained and does not require

% Optimization Toolbox or Statistics and Machine Learning Toolbox.

if isempty(Umat)

    pairs = zeros(0,2);

    return;

end

[K,S] = size(Umat);

if K==0 || S==0

    pairs = zeros(0,2);

    return;

end

% Replace accidental nonfinite utilities by a very unattractive finite value.

U = double(Umat);

finiteVals = U(isfinite(U));

if isempty(finiteVals)

    U(:) = 0;

else

    uMin = min(finiteVals);

    uMax = max(finiteVals);

    penalty = max(1,abs(uMax-uMin)+abs(uMax)+abs(uMin));

    U(~isfinite(U)) = uMin - 1e6*penalty;

end

% The rectangular Hungarian solver below minimizes cost for nRows<=nCols.

% Convert maximum utility to minimum cost.  Orient the smaller side as rows

% so exactly min(K,S) real one-to-one assignments are produced.

maxU = max(U(:));

if K >= S

    % Rows = subchannels, columns = SUs.

    cost = maxU - U.';              % S x K, with S <= K

    colForRow = hungarian_min_rect(cost);

    pairs = zeros(S,2);

    for s = 1:S

        k = colForRow(s);

        pairs(s,:) = [k,s];

    end

else

    % Rows = SUs, columns = subchannels.

    cost = maxU - U;                % K x S, with K < S

    colForRow = hungarian_min_rect(cost);

    pairs = zeros(K,2);

    for k = 1:K

        s = colForRow(k);

        pairs(k,:) = [k,s];

    end

    % Preserve the convention used elsewhere in the code: sorted by SC.

    pairs = sortrows(pairs,2);

end

end
