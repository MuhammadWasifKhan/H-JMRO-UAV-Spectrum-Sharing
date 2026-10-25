function Qnew = resample_single_trajectory_local(Qold,Nnew)

Nold = size(Qold,1);

xOld = linspace(0,1,Nold);

xNew = linspace(0,1,Nnew);

Qnew = zeros(Nnew,3);

Qnew(:,1) = interp1(xOld,Qold(:,1),xNew,'linear','extrap');

Qnew(:,2) = interp1(xOld,Qold(:,2),xNew,'linear','extrap');

Qnew(:,3) = interp1(xOld,Qold(:,3),xNew,'linear','extrap');

end
