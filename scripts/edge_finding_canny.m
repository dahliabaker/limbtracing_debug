function [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid, angdiff)

    %Find edge points
    [ast_edge] = edge(asteroid,'canny',.5,5); 
    %need to flip E_v to match camera orientation (flip vhat 180degrees
    %about uhat axis)
    
    %ast_edge = flip(ast_edge,1);
    E_idx = find(ast_edge);
    
    %[E_u,E_v] = ind2sub( size(asteroid), E_idx);
    [E_v,E_u]  = ind2sub( size(asteroid), E_idx ); 
   
    %calculate middle point (average)
    
%     k = convhull(E_u,E_v);
%     
%     E_u = E_u(k);
%     E_v = E_v(k);

    mid_pt_u = mean(E_u);
    mid_pt_v = mean(E_v);
   
    
    %put all points in frame centered on middle point as origin
    for i = 1:length(E_u)
        centered_u(i) = E_u(i) - mid_pt_u;
        centered_v(i) = E_v(i) - mid_pt_v;
    end
    
    %compute angles between all points and average
    for i = 1:length(centered_u)
       angle(i) = (atan2(centered_v(i),centered_u(i))*180/pi);
    end
    
    for ii = 1:length(angle)
        if (sign(angle(ii)) == -1)
            diff = 180.0+angle(ii);
            angle(ii) = 180+diff;
        end
    end
    points = [centered_u; centered_v];
    B = angle;
    angle_sorted = zeros(size(angle));
    points_sorted = zeros(size(points));
    
    for ii = 1:length(centered_u)
        x = min(B);
        ind = find(B==x);
        if length(ind) > 1
            ind = ind(1);
        end
        B(ind) = NaN;
        angle_sorted(ii) = x;
        points_sorted(1:2,ii) = [centered_u(ind); centered_v(ind)];
    end
    
    j = 1;
    ang = 0;
    for i = 1:length(centered_u)
        if (angle_sorted(i) > ang)
            points_saved(1:2,j) = points_sorted(1:2,i);
            j = j+1;
            ang = ang+angdiff;
        end
    end
    for i = 1:length(points_saved(1,:))
        trim_u(i) = points_saved(1,i);
        trim_v(i) = points_saved(2,i);
    end
    
    
end