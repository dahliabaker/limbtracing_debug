function [intPnt,intDist] = line_patch_intersection(termLine,imgPatch)
    
    % patch normal
    u = (imgPatch(:,1)-imgPatch(:,2));
    v = (imgPatch(:,3)-imgPatch(:,2));
    nPatch = normc(cross(u,v));% unit normal of the patch/plane
    
    % vector in the direction of the line
    uLine = termLine(:,2)-termLine(:,1);
    
    % check if the line is parallel to the patch
    if dot(uLine,nPatch) == 0% check for parallelity
        
        intPnt = [];
        intDist = [];
                
    else% if not parallel, line will intersect the plane containing the patch
        
        % compute intersection distance along the terminator line
        intDist = dot((imgPatch(:,1)-termLine(:,1)),nPatch)/dot(uLine,nPatch);
        
        % compute intersection point along the terminator line
        intPnt = termLine(:,1) + uLine*intDist;
        
        % length of terminator ray
        rayLength = norm(uLine);
        
        % representation of intersection point a basis of vectors in the
        % plane
        intPnt_uv = [dot(intPnt-imgPatch(:,2),normc(u));
                     dot(intPnt-imgPatch(:,2),normc(v))];
        
        % invalidate intersection if any of the checks are tripped
        if intDist > rayLength || intDist < 0 || ...% if point is outside of the terminator line segment
           any(intPnt_uv < 0) || any(intPnt_uv > [norm(u);norm(v)]) % if point is outside the patch boundary
            
            intPnt = [];
            intDist = [];
            
        end
        
    end
    
end