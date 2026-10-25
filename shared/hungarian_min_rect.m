function colForRow = hungarian_min_rect(cost)

% Rectangular Hungarian algorithm for a finite cost matrix with nRows<=nCols.

% Returns colForRow(i)=j, meaning row i is assigned to column j.

cost = double(cost);

[nRows,nCols] = size(cost);

if nRows==0

    colForRow = zeros(0,1);

    return;

end

if nRows > nCols

    error('hungarian_min_rect requires number of rows <= number of columns.');

end

if any(~isfinite(cost(:)))

    error('Hungarian cost matrix must be finite.');

end

% Potentials and augmenting-path arrays. Index 1 plays the role of the

% conventional dummy index 0 in the standard Hungarian formulation.

u   = zeros(nRows+1,1);

v   = zeros(nCols+1,1);

p   = zeros(nCols+1,1);

way = zeros(nCols+1,1);

for i = 1:nRows

    p(1) = i;

    j0 = 1;

    minv = inf(nCols+1,1);

    used = false(nCols+1,1);

    while true

        used(j0) = true;

        i0 = p(j0);

        delta = inf;

        j1 = 0;

        for j = 2:nCols+1

            if ~used(j)

                cur = cost(i0,j-1) - u(i0+1) - v(j);

                if cur < minv(j)

                    minv(j) = cur;

                    way(j) = j0;

                end

                % Deterministic tie breaking by the smallest column index.

                if minv(j) < delta

                    delta = minv(j);

                    j1 = j;

                end

            end

        end

        if ~isfinite(delta) || j1==0

            error('Hungarian assignment failed to find an augmenting path.');

        end

        for j = 1:nCols+1

            if used(j)

                u(p(j)+1) = u(p(j)+1) + delta;

                v(j) = v(j) - delta;

            elseif j>=2

                minv(j) = minv(j) - delta;

            end

        end

        j0 = j1;

        if p(j0)==0

            break;

        end

    end

    % Augment the current matching.

    while true

        j1 = way(j0);

        p(j0) = p(j1);

        j0 = j1;

        if j0==1

            break;

        end

    end

end

colForRow = zeros(nRows,1);

for j = 2:nCols+1

    if p(j)>0

        colForRow(p(j)) = j-1;

    end

end

if any(colForRow==0)

    error('Hungarian assignment returned an incomplete matching.');

end

end

% ========================== KKT power ====================================
