function [trim_u, trim_v] = edge_finding(asteroid,masksize, sample_num, itr_num, angdiff)

    %Read an image
    %asteroid = imread(img);
    %have data about z direction from lbl file

    %Identify masking size % should automate later
    mask = false(size(asteroid));
    
    mask(masksize(1):masksize(2),masksize(3):masksize(4)) = true;

    %Find foreground
    %takes a long time, the 1500 is the minimum iterations required to converge
    bw = activecontour(asteroid, mask, itr_num,'edge');
    %con_idx = find(bw); %con for closed contour
    %[con_v,con_u]  = ind2sub(size(asteroid), con_idx ); 

    %Find edge points
    [ast_edge] = edge(bw,'sobel',.3); 
    E_idx = find(ast_edge);
    
    %[E_u,E_v] = ind2sub( size(asteroid), E_idx);
    [E_v,E_u]  = ind2sub( size(asteroid), E_idx ); 

    %%
    %Sample every x number edge points
%     j = 1;
%     x = sample_num;
%     for i = 1:length(E_u)
%         if mod(i,x) == 0
%             trim_u(j,1) = E_u(i,1);
%             trim_v(j,1) = E_v(i,1);
%             j = j+1;
%         end
%     end
%     trim_u = E_u;
%     trim_v = E_v;
    %sort points in counterclockwise angle order
    
    %calculate middle point (average)
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
        B(ind) = NaN;
        angle_sorted(ii) = x;
        points_sorted(1:2,ii) = [centered_u(ind); centered_v(ind)];
    end
    
    j = 1;
    ang = 0;
    for i = 1:length(centered_u)
        if (angle_sorted(i) > ang)
            %angle_saved(j) = angle_sorted(i);
            points_saved(1:2,j) = points_sorted(1:2,i);
            j = j+1;
            ang = ang+angdiff;
        end
    end
   
    trim_u = points_saved(1,:)';
    trim_v = points_saved(2,:)';
    
    %order from 0 to 360
    %disp('angles created');
    %be done
    
    
    
end