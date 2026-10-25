function d = point_segment_distance_2d(p,a,b)

ab = b-a;

den = dot(ab,ab);

if den <= 1e-16

    d = norm(p-a);

else

    tau = min(max(dot(p-a,ab)/den,0),1);

    d = norm(p-(a+tau*ab));

end

end
