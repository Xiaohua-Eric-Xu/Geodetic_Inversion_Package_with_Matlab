function d = points_to_line(pts, v1, v2)
      a = v2 - v1;
      b = zeros(size(pts));
      d = zeros(size(pts,1),1);
      for i = 1:1:size(pts,1)
        b(i,:) = pts(i,:) - v1;
        d(i) = norm(cross(a,b(i,:))) / norm(a);
      end     
