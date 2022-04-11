function [shapePnts] = shape_from_limb_and_term(limb_starts,limb_ends)
% array lengths
numImages = length(limb_starts);
numLimbRays = length(limb_starts{1})/2;
numTermRays = length(limb_starts{1})/2;

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

%% for each terminator ray, compute intersection point with all patches
shapePntsRaw = cell(numImages,numTermRays);
parfor iImg = 1:numImages
    disp(iImg)
    for iPntT = 1:numTermRays

        % ray-start is closer to the body center than ray-end
        termRay = [limb_ends{iImg}(iPntT+numTermRays,1:3);% ray-start
                   limb_starts{iImg}(iPntT+numTermRays,1:3)]';% ray-end
        
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
                if ~isempty(intDist)
                    
                    % check if the current distance is less than minimum distance
                    if intDist < bestDist
                        bestPnt = intPnt;
                        bestDist = intDist;
                    end
                    
                end
                
            end
        end
        if isempty(bestPnt)
            keyboard
        end
        % save point with the minimum distance as the shape point
        shapePntsRaw{iImg,iPntT} = bestPnt(:)';
        drawnow
    end
    
end
% convert cell array to matrix and remove nans
shapePntsVec = shapePntsRaw(:);
shapePntsMat = cell2mat(shapePntsVec);
shapePnts = shapePntsMat(~isnan(vecnorm(shapePntsMat,2,2)),:);
end