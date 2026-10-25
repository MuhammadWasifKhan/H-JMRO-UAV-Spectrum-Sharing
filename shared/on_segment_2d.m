function tf = on_segment_2d(a,b,p,tol)

tf = p(1)>=min(a(1),b(1))-tol && p(1)<=max(a(1),b(1))+tol && ...
     p(2)>=min(a(2),b(2))-tol && p(2)<=max(a(2),b(2))+tol;

end
