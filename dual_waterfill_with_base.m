function p = dual_waterfill_with_base(a,g,Pmax,Gamma,nu)
% Mode C/C KKT update with individual PU/WU caps and shared-power multiplier.
p = individual_capped_waterfill(a,g,Pmax,Gamma,nu);
end
