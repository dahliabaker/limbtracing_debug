function [shapePnts] = shape_from_limb_and_term(limb_starts,limb_ends,termRays_starts,termRays_ends)
% array lengths
numImages = length(limb_starts);
numLimbRays = length(limb_starts{1});
numTermRays = length(termRays_starts{1});
ax = gca;
%% construct limb patches between pair of neighboring limb rays (along same image)
imgPatch = cell(numImages,numLimbRays-2);
for iImg = 1:numImages
    delete(ax.Children(1:end-144))
    for iPnt = 1:numLimbRays-2
        
        % patch is defined by 4 vertices
        imgPatch{iImg,iPnt} = [limb_starts{iImg}(iPnt,1:3);
                               limb_ends{iImg}(iPnt,1:3);
                               limb_ends{iImg}(iPnt+1,1:3);
                               limb_starts{iImg}(iPnt+1,1:3)]';
        
        patch('XData',imgPatch{iImg,iPnt}(1,:),...
              'YData',imgPatch{iImg,iPnt}(2,:),...
              'ZData',imgPatch{iImg,iPnt}(3,:),...
              'FaceColor','r','FaceAlpha',.1,'LineWidth',1)
        
    end
    
    drawnow
    
end
%% for each terminator ray, compute intersection point with all patches
shapePnts_fromTerm = cell(numImages,numTermRays);
ax = gca;
parfor iImg = 1:numImages
    disp(iImg)
    for iPntT = 1:numTermRays

        % ray-start is closer to the body center than ray-end
        termRay = [termRays_ends{iImg}(iPntT,1:3);% ray-start
                   termRays_starts{iImg}(iPntT,1:3)]';% ray-end
%         if exist('tRay','var')
%             delete(tRay)
%         end
%         tRay = plot3(termRay(1,:),termRay(2,:),termRay(3,:),'Color','k','Marker','*','LineStyle','-.','LineWidth',2);
%         drawnow
        % loop through limb patches and compute intersection with current terminator ray
        bestPnt = nan(3,1); bestDist = Inf;bestPatch_idx = [1,1];
        for ii = 1:size(imgPatch,1)
            
            % skip iteration if limb image set is same as terminator image set
            if ii == iImg
                continue
            end
            for jj = 1:size(imgPatch,2)
                
                % compute intersection between terminator ray and limb patch
                [intPnt,intDist] = line_patch_intersection(termRay,imgPatch{ii,jj});
                if ~isempty(intDist)
                    % check if the current distance is less than minimum distance
                    if intDist < bestDist
                        bestPnt = intPnt;
                        bestDist = intDist;
                        bestPatch_idx = [ii,jj];
                    end
                    
                end
                
            end
        end
%         delete(ax.Children(1:end-144))
%         plot3(bestPnt(1),bestPnt(2),bestPnt(3),'sg','MarkerFaceColor','g');
%         patch('XData',imgPatch{bestPatch_idx(1),bestPatch_idx(2)}(1,:),...
%               'YData',imgPatch{bestPatch_idx(1),bestPatch_idx(2)}(2,:),...
%               'ZData',imgPatch{bestPatch_idx(1),bestPatch_idx(2)}(3,:),...
%               'FaceColor','r','FaceAlpha',.1,'LineWidth',1),drawnow
          
        % save point with the minimum distance as the shape point
        shapePnts_fromTerm{iImg,iPntT} = bestPnt(:)';
        
    end
    
end
% convert cell array to matrix and remove nans
shapePntsVec = shapePnts_fromTerm(:);
shapePntsMat = cell2mat(shapePntsVec);
shapePnts = shapePntsMat(~isnan(vecnorm(shapePntsMat,2,2)),:);
end