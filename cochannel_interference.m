function I = cochannel_interference(idxList,p,hProtect)
% I(ss) is interference at the receiver assigned to subchannel ss.
I=zeros(numel(hProtect),1);
for ell=1:size(idxList,1)
    ss=idxList(ell,2);
    I(ss)=I(ss)+max(p(ell),0)*max(hProtect(ss),0);
end
end
