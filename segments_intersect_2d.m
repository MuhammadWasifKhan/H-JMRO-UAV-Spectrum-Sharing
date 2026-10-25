function tf = segments_intersect_2d(a,b,c,d)

tol = 1e-10;

o1 = orient2d(a,b,c); o2 = orient2d(a,b,d);

o3 = orient2d(c,d,a); o4 = orient2d(c,d,b);

tf = ((o1>tol && o2<-tol) || (o1<-tol && o2>tol)) && ...
     ((o3>tol && o4<-tol) || (o3<-tol && o4>tol));

if tf, return; end

if abs(o1)<=tol && on_segment_2d(a,b,c,tol), tf=true; return; end

if abs(o2)<=tol && on_segment_2d(a,b,d,tol), tf=true; return; end

if abs(o3)<=tol && on_segment_2d(c,d,a,tol), tf=true; return; end

if abs(o4)<=tol && on_segment_2d(c,d,b,tol), tf=true; return; end

tf=false;

end
