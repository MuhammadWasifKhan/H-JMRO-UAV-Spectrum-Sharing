function [Qprop,convHistory,Qiter] = run_joint_traj_quiet_3uav(varargin)
% Suppress diagnostic output without changing any AO/SCA computations.
evalc('[Qprop,convHistory,Qiter] = build_joint_traj_AO_actual_3uav(varargin{:});');
end
