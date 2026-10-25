function [Qcand,statusNow] = solve_traj_cvx_subproblem_3uav( ...
        Qlocal,scaData,q_start_all,q_end_all,N,U, ...
        localTrust,stepMax,marginFrac,GammaLic,GammaUn,area)

% =========================================================================
% JOURNAL-STYLE JOINT CVX TRAJECTORY SUBPROBLEM
% =========================================================================
% Same logic used by the working 2-C-UAV implementation:
%   maximize scaled first-order rate improvement + small proximal term,
%   subject only to trajectory geometry/mobility/separation/trust constraints.
% The normal communication optimizer is re-solved after CVX and re-enforces
% individual PU/WU interference protection at the candidate trajectory.
%
% marginFrac/GammaLic/GammaUn remain in the signature only so the surrounding
% ICC code interface is unchanged; they are intentionally not used here.
% =========================================================================

global UAV_HMIN UAV_HMAX UAV_DMIN UAV_DMIN_XY UAV_XMIN UAV_XMAX

Qcand = Qlocal;
statusNow = 'Failed';

if exist('cvx_begin','file') ~= 2
    statusNow = 'CVX not installed';
    return;
end

% Build one N-by-(3U) gradient matrix from the already optimized/frozen
% mode, Hungarian and KKT solution at the current expansion point.
X0 = zeros(N,3*U);
G = zeros(N,3*U);

for n = 1:N
    for uu = 1:U
        cols = (3*uu-2):(3*uu);
        X0(n,cols) = Qlocal(n,:,uu);
    end

    DL = scaData{n}.lic;
    for rr = 1:numel(DL)
        uu = DL(rr).u;
        cols = (3*uu-2):(3*uu);
        G(n,cols) = G(n,cols) + DL(rr).gradRate(:).';
    end

    DU = scaData{n}.un;
    for rr = 1:numel(DU)
        uu = DU(rr).u;
        cols = (3*uu-2):(3*uu);
        G(n,cols) = G(n,cols) + DU(rr).gradRate(:).';
    end
end

% Numerical scaling exactly follows the working journal-style logic.
scale = max(1,max(abs(G(:))));

% Explicitly clear a CVX model left in THIS function workspace by a prior
% failed/interrupted CVX construction. cvx_clear alone may not remove the
% local cvx_problem variable that triggers the overwrite warning.
clear cvx_problem
if exist('cvx_clear','file') == 2
    cvx_clear;
end

try
    cvx_begin quiet
        variable X(N,3*U)

        maximize( sum(sum((G/scale).*(X-X0))) - ...
            (1e-4/max(localTrust^2,1e-9))* ...
            sum_square(reshape(X-X0,N*3*U,1)) )

        subject to
            % Endpoints, region/altitude bounds, trust region and mobility.
            for uu = 1:U
                cols = (3*uu-2):(3*uu);

                X(1,cols) == q_start_all(uu,:);
                X(N,cols) == q_end_all(uu,:);

                for n = 1:N
                    X(n,cols(1)) >= UAV_XMIN(uu);
                    X(n,cols(1)) <= UAV_XMAX(uu);
                    X(n,cols(2)) >= 0;
                    X(n,cols(2)) <= area;
                    X(n,cols(3)) >= UAV_HMIN;
                    X(n,cols(3)) <= UAV_HMAX;

                    norm(X(n,cols)-X0(n,cols),2) <= localTrust;

                end

                for n = 1:N-1
                    norm(X(n+1,cols)-X(n,cols),2) <= stepMax;
                end
            end

            % First-order individual PU/WU cochannel interference constraints.
            % Every allocated subchannel protects its corresponding receiver.
            % Nonlinear feasibility is checked after resource reoptimization.
            for n = 1:N
                PL = scaData{n}.protLic;
                for rr = 1:numel(PL)
                    uu = PL(rr).u;
                    cols = (3*uu-2):(3*uu);
                    PL(rr).I0 + PL(rr).gradI(:).' * ...
                        (X(n,cols)-X0(n,cols)).' <= GammaLic;
                end
                PU = scaData{n}.protUn;
                for rr = 1:numel(PU)
                    uu = PU(rr).u;
                    cols = (3*uu-2):(3*uu);
                    PU(rr).I0 + PU(rr).gradI(:).' * ...
                        (X(n,cols)-X0(n,cols)).' <= GammaUn;
                end
            end

            % Conservative first-order inner approximation of synchronized
            % 3-D and horizontal separation constraints.
            for n = 1:N
                for uu = 1:U-1
                    colsU = (3*uu-2):(3*uu);
                    for vv = uu+1:U
                        colsV = (3*vv-2):(3*vv);

                        d0 = X0(n,colsU)-X0(n,colsV);
                        if sum(d0.^2) < UAV_DMIN^2-1e-5
                            error(['Journal-style SCA requires a feasible ' ...
                                   '3-D separation expansion point.']);
                        end
                        2*d0*(X(n,colsU)-X(n,colsV)).' - ...
                            sum(d0.^2) >= UAV_DMIN^2;

                        d0xy = X0(n,colsU(1:2))-X0(n,colsV(1:2));
                        if sum(d0xy.^2) < UAV_DMIN_XY^2-1e-5
                            error(['Journal-style SCA requires a feasible ' ...
                                   'horizontal separation expansion point.']);
                        end
                        2*d0xy*(X(n,colsU(1:2))-X(n,colsV(1:2))).' - ...
                            sum(d0xy.^2) >= UAV_DMIN_XY^2;
                    end
                end
            end
    cvx_end

    statusNow = cvx_status;

    if strcmpi(cvx_status,'Solved') || strcmpi(cvx_status,'Inaccurate/Solved')
        for uu = 1:U
            cols = (3*uu-2):(3*uu);
            Qcand(:,:,uu) = full(X(:,cols));
        end
    else
        Qcand = Qlocal;
    end
    clear cvx_problem

catch ME
    statusNow = ['Failed: ',ME.message];
    Qcand = Qlocal;
    clear cvx_problem
    if exist('cvx_clear','file') == 2
        try
            cvx_clear;
        catch
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
