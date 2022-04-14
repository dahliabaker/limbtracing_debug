function [shapeEndPnts, shapePnts, shapePntNhats] = shape_from_limbs_lo(limbRayStarts, limbRayEnds, longitudeSet, numPnts)
%

%if limb_only = 0, then it's limb and terminator
%if limb_only = 1, then it's only limbs
% Remember limbs are not necessarily planar! 
%
% If connect two neighboring points from a limb, then make it a plane in
% the viewing direction, can get intersection of planes, and those lines
% are lines on the shape model
%   -> this technically requires my limbs to have uhat info, but to start I
%   am assuming that all "pixels" are viewing in the identical uhat

% Test this with a couple (2 then 3) circle, triangle or square limbs

% NOT JUST INTERSECTIONS OF LIMB PATCHES ON SHAPE, BUT EVERYTHING BETWEEN
% THEM!
%  -> I'll have to keep track of the intersection lines on each patch and
%  then populate the area between them...


% Angle threshold for being "on" a line
angThresh = 1e-2; % deg
cangThresh = cosd(angThresh);

numLimbPatches = size(limbRayStarts{1},1); 
totalLimbPatches = numLimbPatches*length(longitudeSet); 

% have to save totalLimbPatches of each output
% Indices: organized in groups by longitude, with the first one having
% numLimbPatches-1 entries and each subsequent group having one less
% Its the super-diagonal portion of an totalLimbPatches x totalLimbPatches
% matrix
% Just doing the 2D way for now..   

shapeEndPnts = cell(totalLimbPatches,totalLimbPatches); %(totalLimbPatches*(totalLimbPatches-1))/2,1);
% shape_tInt = cell(totalLimbPatches,totalLimbPatches);
% shape_uhat = cell(totalLimbPatches,totalLimbPatches);
% shape_x0 = cell(totalLimbPatches,totalLimbPatches);
shape_nhat = cell(totalLimbPatches,1); % different, just the normal for each patch

numEndPnts = 2; % this will just return the endpoints, which I think is all i need right now

skipmm = [];

% Compute all the intersections between limb patches
for ii = 1:length(longitudeSet)-1
    disp(ii)
    for jj = 1:numLimbPatches 
        
        if jj < numLimbPatches
            limb1 = [limbRayStarts{ii}(jj,:); limbRayStarts{ii}(jj+1,:); limbRayEnds{ii}(jj,:); limbRayEnds{ii}(jj+1,:)]'; 
%         elseif limb_only == 0
%             limb1 = [limbRayStarts{ii}(jj,:); limbRayStarts{ii}(1,:); limbRayEnds{ii}(jj,:); limbRayEnds{ii}(1,:)]'; 
%         
        end
        
        for mm = ii+1:length(longitudeSet)
            
            %skip if limbs are 180 degrees apart
%             if abs(longitudeSet(mm) - longitudeSet(ii) - pi) < pi/180
%                 if ii == 1
%                     skipmm = [skipmm; mm];
%                 end
%                 continue
%             end
            
            for nn = 1:numLimbPatches 
                
                if nn < numLimbPatches
                    limb2 = [limbRayStarts{mm}(nn,:); limbRayStarts{mm}(nn+1,:); limbRayEnds{mm}(nn,:); limbRayEnds{mm}(nn+1,:)]';
%                 elseif limb_only == 0
%                     limb2 = [limbRayStarts{mm}(nn,:); limbRayStarts{mm}(1,:); limbRayEnds{mm}(nn,:); limbRayEnds{mm}(1,:)]';
                end
                
%                 [shapeEndPnts{(ii-1)*numLimbPatches + jj,(mm-1)*numLimbPatches + nn}, shape_tInt{(ii-1)*numLimbPatches + jj,(mm-1)*numLimbPatches + nn}, shape_uhat{(ii-1)*numLimbPatches + jj,(mm-1)*numLimbPatches + nn}, shape_x0{(ii-1)*numLimbPatches + jj,(mm-1)*numLimbPatches + nn}, nhat1, nhat2] = find_shape_pts_from_limb_segments(limb1, limb2, numEndPnts);
                [shapeEndPnts{(ii-1)*numLimbPatches + jj,(mm-1)*numLimbPatches + nn}, ~, ~, ~, nhat1, nhat2] = find_shape_pts_from_limb_segments(limb1, limb2, numEndPnts);
                if ii == 1
                    shape_nhat{(ii-1)*numLimbPatches + jj} = nhat1;
                    if jj == 1
                        shape_nhat{(mm-1)*numLimbPatches + nn} = nhat2;
                    end
                elseif ~isempty(find(skipmm==mm,1))
                    if jj == 1
                        shape_nhat{(mm-1)*numLimbPatches + nn} = nhat2;
                        if nn == numLimbPatches
                            skipmm(skipmm==mm) = [];
                        end
                    end
                elseif ~isempty(find(skipmm==ii,1))
                    shape_nhat{(ii-1)*numLimbPatches + jj} = nhat1;
                    if jj == numLimbPatches
                        skipmm(skipmm==ii) = [];
                    end
                end
                
            end
            
        end
        
    end
    
end

% On each patch, want to whittle away, keeping the smallest segment of that
% patch
% To determine which one is "inside", project the inward normal of the
% intersecting patch onto the patch being whittled. If two neighboring
% intersection points have these projected normals facing toward each
% other, they are opposite ends and should both be kept. If the neighboring
% points have the normals in the same direction, then the one that is
% further outside with respect to this normal direction should be dropped. 

% At first, maybe only keep intersection points on the edges of the patch?
% This would decrease accuracy, but be simpler to code up
% Yes, start this way - only worry about cutting down each ray

% Possibly easier to process full longitude sets together? or is fully
% sequential fine such that order doesn't matter?
% HAVE TO BATCH ALL INTERSECTIONS OF A RAY in order to catch concavities


% HOW TO FIX ITOKAWA:
% Instead of adding portions to the array as I did before, for each trimming
% view (other views) I need to compute the segment that is inside that
% viewing cone and then take the intersection with what I already had
%
% Functionally this means that I construct a rayIntersections for each view
% (ii of longitudeSet), find the portions of the ray to be included from
% that view, and then take the intersection of those sets with the overall
% set
%%%%%%%Dahlia's Ray Length Statistics Edit%%%%%%%%%%%%%%%%%%%%%
%first calculate the stats of the rays
% rayLengths = [];
% % 
shapePnts = zeros(0,3);
shapePntNhats = zeros(0,3);
for ii = 1:length(longitudeSet)
    disp(ii)
    for jj = 1:numLimbPatches
        
%         if (ii==2 && jj == 23) || ii==4 || ii==10 %ii == 12 && (jj == 5) % 11, 5 itokawa
%             disp('hi')
%         end
        
        start2end = (limbRayEnds{ii}(jj,:) - limbRayStarts{ii}(jj,:))';
        rayLength = norm(start2end);
        rayInd = (ii-1)*numLimbPatches + jj;
        
        
        
        keepSegments = [0 1]; % c1 = %start, c2 = %stop; add rows if new sections;
        
        otherViews = 1:length(longitudeSet);
        otherViews(ii) = [];
        for kk = otherViews
            
            % For each view, find segments which are inside this view's
            % viewing cone
            % Column 1 is the percent in the ray
            % Column 2 is the normal direction compared to the ray direction
            %  -> start with end points of the ray 
            rayIntersections = [0 -1; 1 1];
            
            for mm = 1:numLimbPatches
        
                % Find where this patch trims the ray and its projected
                % normal
                pnum = (kk-1)*numLimbPatches + mm;
                if (rayInd < pnum && ~isempty(shapeEndPnts{rayInd,pnum})) || (rayInd > pnum && ~isempty(shapeEndPnts{pnum,rayInd}))
                    for nn = 1:2
                        if rayInd < pnum
                            pt = shapeEndPnts{rayInd,pnum}(:,nn);
                        else
                            pt = shapeEndPnts{pnum,rayInd}(:,nn);
                        end
                        start2pt = pt - limbRayStarts{ii}(jj,:)';
                        s2pDist = norm(start2pt);
                        
                        % check to see if this point is on this ray
                        angArg = dot(start2pt,start2end)/rayLength/s2pDist;
                        
                        if angArg > cangThresh
                            % get percentage of ray
                            pcntIn = s2pDist/rayLength;
                            
                            % compare percentage and normal to determine if should be
                            % inserted into array
                            projNormal = dot(start2end./rayLength, shape_nhat{pnum}(:,1));
                            
                            rayIntersections = [rayIntersections; pcntIn sign(projNormal)];
                            
                        end
                        
                    end
                    
                end
                
            end
            
            % Determine the segment(s) of the ray within this viewing cone
            % Sort through the intersections and keep appropriate sets
        
            % Sort by pcntIn column
            rayIntersections = sortrows(rayIntersections);
            
            % Throw out any 0 entries for normal directions if they exist
            zeroNormInd = rayIntersections(:,2)==0;
            rayIntersections(zeroNormInd,:) = [];
            
            % Find change in signs
            normalChange = diff(rayIntersections(:,2));
            includeInds = find(normalChange==2);
            
            % set of current segments
            currSegs = zeros(length(includeInds),2);
            for qq = 1:length(includeInds)
                currSegs(qq,:) = [rayIntersections(includeInds(qq),1), rayIntersections(includeInds(qq)+1,1)];
            end
            
            % take the intersection of this set of segments with the
            % running set
            intersectSegments = zeros(0,2);
            for mm = 1:size(keepSegments,1)
                for nn = 1:size(currSegs,1)
                    % Keep the intersection of the arcs that intersects both patches
                    t_intersect = [max([keepSegments(mm,1) currSegs(nn,1)]) min([keepSegments(mm,2) currSegs(nn,2)])];
                    
                    if t_intersect(1) <= t_intersect(2)
                        intersectSegments = [intersectSegments; t_intersect(1), t_intersect(2)];
                    end
                end
            end
            
            keepSegments = intersectSegments;
            
        end
        
        
        if isempty(keepSegments) || (keepSegments(1,1) == 0 && keepSegments(1,2) == 1)
            %fprintf('Threw out untrimmed ray index %d (long %d, az %d)\n',rayInd,ii,jj)
            keepSegments = zeros(0,2);
        end

            
        for kk = 1:size(keepSegments,1)
            pntPcnts = linspace(keepSegments(kk,1),keepSegments(kk,2),numPnts);
            pnts = limbRayStarts{ii}(jj,:)' + start2end.*pntPcnts;
            shapePnts = [shapePnts; pnts'];
            nh1 = shape_nhat{rayInd}(:,1);
            if jj == 1
                nh2 = shape_nhat{ii*numLimbPatches}(:,1);
            else
                nh2 = shape_nhat{rayInd-1}(:,1);
            end
            nhat = (nh1'+nh2')./2;
            nhat = nhat./norm(nhat);
            shapePntNhats = [shapePntNhats; repmat(nhat,numPnts,1)];
        end       
    end
end

% figure
% histogram(rayLengths,100)
% title('Post-Trim Ray Length Distribution','FontSize', 18)
% xlabel('Ray Length (km)','FontSize', 16)
% ylabel('Count','FontSize', 16)
% for ii = 1:length(longitudeSet)
%     for jj = 1:numLimbPatches
%         start2end = (limbRayEnds{ii}(jj,:) - limbRayStarts{ii}(jj,:))';
%         rayLength = norm(start2end);
%         rayInd = (ii-1)*numLimbPatches + jj;
% 
% 
% 
%         keepSegments = [0 1]; % c1 = %start, c2 = %stop; add rows if new sections;
% 
%         otherViews = 1:length(longitudeSet);
%         otherViews(ii) = [];
%         for kk = otherViews
% 
%             % For each view, find segments which are inside this view's
%             % viewing cone
%             % Column 1 is the percent in the ray
%             % Column 2 is the normal direction compared to the ray direction
%             %  -> start with end points of the ray 
%             rayIntersections = [0 -1; 1 1];
% 
%             for mm = 1:numLimbPatches
% 
%                 % Find where this patch trims the ray and its projected
%                 % normal
%                 pnum = (kk-1)*numLimbPatches + mm;
%                 if (rayInd < pnum && ~isempty(shapeEndPnts{rayInd,pnum})) || (rayInd > pnum && ~isempty(shapeEndPnts{pnum,rayInd}))
%                     for nn = 1:2
%                         if rayInd < pnum
%                             pt = shapeEndPnts{rayInd,pnum}(:,nn);
%                         else
%                             pt = shapeEndPnts{pnum,rayInd}(:,nn);
%                         end
%                         start2pt = pt - limbRayStarts{ii}(jj,:)';
%                         s2pDist = norm(start2pt);
% 
%                         % check to see if this point is on this ray
%                         angArg = dot(start2pt,start2end)/rayLength/s2pDist;
% 
%                         if angArg > cangThresh
%                             % get percentage of ray
%                             pcntIn = s2pDist/rayLength;
% 
%                             % compare percentage and normal to determine if should be
%                             % inserted into array
%                             projNormal = dot(start2end./rayLength, shape_nhat{pnum}(:,1));
% 
%                             rayIntersections = [rayIntersections; pcntIn sign(projNormal)];
% 
%                         end
% 
%                     end
% 
%                 end
% 
%             end
% 
%             % Determine the segment(s) of the ray within this viewing cone
%             % Sort through the intersections and keep appropriate sets
% 
%             % Sort by pcntIn column
%             rayIntersections = sortrows(rayIntersections);
% 
%             % Throw out any 0 entries for normal directions if they exist
%             zeroNormInd = rayIntersections(:,2)==0;
%             rayIntersections(zeroNormInd,:) = [];
% 
%             % Find change in signs
%             normalChange = diff(rayIntersections(:,2));
%             includeInds = find(normalChange==2);
% 
%             % set of current segments
%             currSegs = zeros(length(includeInds),2);
%             for qq = 1:length(includeInds)
%                 currSegs(qq,:) = [rayIntersections(includeInds(qq),1), rayIntersections(includeInds(qq)+1,1)];
%             end
% 
%             % take the intersection of this set of segments with the
%             % running set
%             intersectSegments = zeros(0,2);
%             for mm = 1:size(keepSegments,1)
%                 for nn = 1:size(currSegs,1)
%                     % Keep the intersection of the arcs that intersects both patches
%                     t_intersect = [max([keepSegments(mm,1) currSegs(nn,1)]) min([keepSegments(mm,2) currSegs(nn,2)])];
% 
%                     if t_intersect(1) <= t_intersect(2)
%                         intersectSegments = [intersectSegments; t_intersect(1), t_intersect(2)];
%                     end
%                 end
%             end
% 
%             keepSegments = intersectSegments;
% 
%         end
% 
% 
%         if isempty(keepSegments) || (keepSegments(1,1) == 0 && keepSegments(1,2) == 1)
%             %fprintf('Threw out untrimmed ray index %d (long %d, az %d)\n',rayInd,ii,jj)
%             keepSegments = zeros(0,2);
%         end
% 
% 
%         for kk = 1:size(keepSegments,1)
%             pntPcnts = linspace(keepSegments(kk,1),keepSegments(kk,2),numPnts);
%             pnts = limbRayStarts{ii}(jj,:)' + start2end.*pntPcnts;
%             %if norm(pnts(:,end)-pnts(:,1)) < (mean(rayLengths)+stdnum*std(rayLengths))
% 
%                 
%                 nh1 = shape_nhat{rayInd}(:,1);
%                 if jj == 1
%                     nh2 = shape_nhat{ii*numLimbPatches}(:,1);
%                 else
%                     nh2 = shape_nhat{rayInd-1}(:,1);
%                 end
%                 nhat = (nh1'+nh2')./2;
%                 nhat = nhat./norm(nhat);
%                 %vec_diff = nhat-spin_pole;
% %                 shapePnts = [shapePnts; pnts'];
% %                 shapePntNhats = [shapePntNhats; repmat(nhat,numPnts,1)];
% %                 if vec_diff(3) < 1.1e-2 || vec_diff(3) > (2-1.1e-2)
% %                     skip this ray, too close to spin pole
% %                     remove these two lines later, just taking out the pole
% %                     adjustment functionality to test
% %                     shapePnts = [shapePnts; pnts'];
% %                     shapePntNhats = [shapePntNhats; repmat(nhat,numPnts,1)];
% %                 else
%                 shapePnts = [shapePnts; pnts'];
%                 shapePntNhats = [shapePntNhats; repmat(nhat,numPnts,1)];
%                 %end
%             %end
%         end
%     end
% end
% 
% 
% 
% 
% % for ii = 1:length(longitudeSet)
% %     for jj = 1:numLimbPatches
% %         intersectSegments = keepSegments{ii,jj};
% %         % With this discretization scheme, and only getting limbs from the
% %         % equator, you will pretty much always have one ray at each pole
% %         % that isn't intersected; we will throw this ray out
% %         
% %         
% %         % Compute the points to keep from these segments
% %         for kk = 1:size(intersectSegments,1)
% %             
% %             pntPcnts = linspace(intersectSegments(kk,1),intersectSegments(kk,2),numPnts);
% %             pnts = limbRayStarts{ii}(jj,:)' + start2end.*pntPcnts;
% %             if norm(pnts(:,end)-pnts(:,1)) < (mean(rayLengths)+2*std(rayLengths))
% %             
% %                 shapePnts = [shapePnts; pnts'];
% %                 nh1 = shape_nhat{rayInd}(:,1);
% %                 if jj == 1
% %                     nh2 = shape_nhat{ii*numLimbPatches}(:,1);
% %                 else
% %                     nh2 = shape_nhat{rayInd-1}(:,1);
% %                 end
% %                 nhat = (nh1'+nh2')./2;
% %                 nhat = nhat./norm(nhat);
% %                 shapePntNhats = [shapePntNhats; repmat(nhat,numPnts,1)];
% %             end
% %             
% %         end
% %     end
% % end


% remove points that are copies of one another
[shapePnts, IA, IC] = uniquetol(shapePnts,'ByRows',true);
shapePntNhats = shapePntNhats(IA,:);
end
% add points in limb plane patches between the ray points computed above


% save data to file
% fid = fopen([outFileName '.xyz'],'w');
% for ii = 1:length(shapePnts)
%     fprintf(fid,'%f %f %f %f %f %f\n',shapePnts(ii,1),shapePnts(ii,2),shapePnts(ii,3), shapePntNhats(ii,1), shapePntNhats(ii,2), shapePntNhats(ii,3));
% end
% fclose(fid)
                