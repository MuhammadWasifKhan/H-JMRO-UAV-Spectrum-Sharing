function d = segment_segment_distance_2d(p1,p2,q1,q2)

% Exact Euclidean distance between two closed line segments in R^2.

if segments_intersect_2d(p1,p2,q1,q2)

    d = 0;

    return;

end

d = min([point_segment_distance_2d(p1,q1,q2), ...
         point_segment_distance_2d(p2,q1,q2), ...
         point_segment_distance_2d(q1,p1,p2), ...
         point_segment_distance_2d(q2,p1,p2)]);

end
