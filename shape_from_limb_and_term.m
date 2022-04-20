function [shapePnts] = shape_from_limb_and_term(limb_starts,limb_ends,termRays_starts,termRays_ends)
% array lengths
numImages = length(limb_starts);
numLimbRays = length(limb_starts{1});
numTermRays = length(termRays_starts{1});
ax = gca;
%% construct limb patches between pair of neighboring limb rays (along same image)
imgPatch = cell(numImages,numLimbRays-1);
for iImg = 1:numImages
    
    for iPnt = 1:numLimbRays-1
        
        % patch is defined by 4 vertices
        imgPatch{iImg,iPnt} = [limb_starts{iImg}(iPnt,1:3);
                               limb_ends{iImg}(iPnt,1:3);
                               limb_ends{iImg}(iPnt+1,1:3);
                               limb_starts{iImg}(iPnt+1,1:3)]';
    end
    
end

% pretrim patches against each other:
for ii = 1:size(imgPatch,1)
    for jj = 1:size(imgPatch,2)
        patch1 = imgPatch{ii,jj};
%         patch('XData',patch1(1,:),...
%               'YData',patch1(2,:),...
%               'ZData',patch1(3,:),...
%               'FaceColor','b','FaceAlpha',.7,'LineWidth',1)
          
        % only trim current patch against immediate neighbors
%         if ii == 1
%             kk = [size(imgPatch,1),2];
%         elseif ii == size(imgPatch,1)
%             kk = [1,size(imgPatch,1)-1];
%         else
%             kk = ii+[-1,1];
%         end
        kk = 1:numImages;kk(ii) = [];
        
        intPatch = [];
        for k = kk
            for nn = 1:size(imgPatch,2)
                patch2 = imgPatch{k,nn};
%                 patch('XData',patch2(1,:),...
%                       'YData',patch2(2,:),...
%                       'ZData',patch2(3,:),...
%                       'FaceColor','r','FaceAlpha',.1,'LineWidth',1)
                [intLines_temp] = find_shape_pts_from_limb_segments(patch1(:,[1,2,4,3]),...
                                                                            patch2(:,[1,2,4,3]),...
                                                                            2);
                if ~isempty(intLines_temp)

%                     plot3(intLines_temp(1,:),intLines_temp(2,:),intLines_temp(3,:),'*g')
                end
                intPatch = [intPatch,intLines_temp];
            end
        end
        if ~isnan(norm(intPatch))
            intPatchU = uniquetol(intPatch','ByRows',true)';
%             plot(alphaShape(intPatch'));
            patch('XData',intPatchU(1,:),...
                  'YData',intPatchU(2,:),...
                  'ZData',intPatchU(3,:),...
                  'FaceColor','r','FaceAlpha',.7,'LineWidth',1)
        end
    end
end

%% for each terminator ray, compute intersection point with all patches
shapePnts_fromTerm = cell(numImages,numTermRays);
% ax = gca;
for iImg = 50:numImages
    disp(iImg)
    for iPntT = 1:numTermRays

        % ray-start is closer to the body center than ray-end
        termRay = [termRays_ends{iImg}(iPntT,1:3);% ray-start
                   termRays_starts{iImg}(iPntT,1:3)]';% ray-end
%         if exist('tRay','var')
%             delete(tRay)
%         end
%         tRay = plot3(termRay(1,:),termRay(2,:),termRay(3,:),'Color','k','Marker','*','LineStyle','-.','LineWidth',2);
        drawnow
        % loop through limb patches and compute intersection with current terminator ray
        bestPnt = nan(3,1); bestDist = Inf;
        for ii = 1:size(imgPatch,1)
            
            % skip iteration if limb image set is same as terminator image set
            if ii == iImg
                continue
            end
            for jj = 1:size(imgPatch,2)
                
                % compute intersection between terminator ray and limb patch
                [intPnt,intDist] = line_patch_intersection(termRay,imgPatch{ii,jj});
%                 p(idx) = patch('XData',imgPatch{ii,jj}(1,:),...
%                     'YData',imgPatch{ii,jj}(2,:),...
%                     'ZData',imgPatch{ii,jj}(3,:),...
%                     'FaceColor','r','FaceAlpha',.1,'LineWidth',1);idx = idx+1;drawnow
                if ~isempty(intDist)
                    % check if the current distance is less than minimum distance
                    if intDist < bestDist
                        bestPnt = intPnt;
                        bestDist = intDist;
%                         if exist('bestP','var')
%                             delete(bestP)
%                         end
                        bestP = plot3(bestPnt(1),bestPnt(2),bestPnt(3),'sr','MarkerFaceColor','r');drawnow
                    end
                    
                end
                
            end
        end
%         plot3(bestPnt(1),bestPnt(2),bestPnt(3),'sr','MarkerFaceColor','r');drawnow
        % save point with the minimum distance as the shape point
        shapePnts_fromTerm{iImg,iPntT} = bestPnt(:)';
        
    end
    
end
% convert cell array to matrix and remove nans
shapePntsVec = shapePnts_fromTerm(:);
shapePntsMat = cell2mat(shapePntsVec);
shapePnts = shapePntsMat(~isnan(vecnorm(shapePntsMat,2,2)),:);
end