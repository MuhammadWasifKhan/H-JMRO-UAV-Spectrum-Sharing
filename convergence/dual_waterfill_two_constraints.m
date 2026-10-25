function p = dual_waterfill_two_constraints(a,g,Pmax,Gamma)
% Compatibility wrapper: individual cochannel protection, not aggregate.
p = individual_capped_waterfill(a,g,Pmax,Gamma,0);
end
