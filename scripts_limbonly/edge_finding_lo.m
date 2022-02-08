function [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_lo(asteroid)

    %Find edge points
    [ast_edge] = edge(asteroid,'canny',.3,3); 
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

        %alphaShape edits
    
    alp = alphaShape(E_v,E_u);
    E_v = alp.Points(:,1);
    E_u = alp.Points(:,2);


%     mid_pt_u = mean(E_u);
%     mid_pt_v = mean(E_v);
    mid_pt_u = 512;
    mid_pt_v = 512;
    
    
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
        
%used to be points_saved, switched to points_sorted 10/8
    for i = 1:length(points_sorted(1,:))
        trim_u(i) = points_sorted(1,i);
        trim_v(i) = points_sorted(2,i);
    end
    
    
end