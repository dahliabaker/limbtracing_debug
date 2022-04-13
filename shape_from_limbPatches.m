% function shapePnts = shape_from_limbPatches(varargin)
% 
% if nargin==1% we are passing in the limb patches
%     imgPatch = varargin{1};% will it be a cell array?
%     
%     % local variables
%     numImages = size(imgPatch,1);
%     numLimbRays = size(imgPatch,2)+1;
%     
% elseif nargin == 2% we are passing in the raw limbs and need to construct the patches
    
    % construct the limb patches
%     limb_starts = varargin{1};
%     limb_ends   = varargin{2};
    
    % local variables
    numImages = length(limb_starts);
    numLimbRays = length(limb_starts{1});
    
    % initialize
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
% else
%     error('incorrect number of inputs')
% end
    
shapePntsRaw = cell(size(imgPatch));
for ii = 1:size(imgPatch,1)
    disp(ii)
    for jj = 1:size(imgPatch,2)
        
        patch1 = imgPatch{ii,jj};
%         patch('XData',patch1(1,:),'YData',patch1(2,:),'ZData',patch1(3,:),...
%             'FaceColor','r','FaceAlpha',.2,'LineWidth',1)
        
        intLines = cell(size(imgPatch));
        shapePnts_temp = [];
        for kk = 1:size(imgPatch,1)
            if kk == ii
                continue% skip intersection between same image sets
            end
            for nn = 1:size(imgPatch,2)
                patch2 = imgPatch{kk,nn};
                [intLines_temp] = find_shape_pts_from_limb_segments(patch1(:,[1,2,4,3]),...
                                                                    patch2(:,[1,2,4,3]),...
                                                                    2);
%                 patch('XData',patch2(1,:),'YData',patch2(2,:),'ZData',patch2(3,:),...
%                     'FaceColor','b','FaceAlpha',.2,'LineWidth',1)
                
                
                % append shape pts
                if ~isempty(intLines_temp)
                    shapePnts_temp = [shapePnts_temp,intLines_temp];
                    intLines{kk,nn} = intLines_temp;
                                       plot3(intLines_temp(1,:),intLines_temp(2,:),intLines_temp(3,:),...
                                             'sr','MarkerFaceColor','r')
                    %                    keyboard
                end
            end
        end
        % find point that is closest to the center
        [minDist,iMin] = min(vecnorm(shapePnts_temp));
        shapePntsRaw{ii,jj} = shapePnts_temp(:,iMin);
        plot3(shapePnts_temp(1,iMin),shapePnts_temp(2,iMin),shapePnts_temp(3,iMin),'sb','MarkerFaceColor','b');drawnow
        keyboard
    end
    
    
end


%% convert cell array to matrix and remove nans
shapePntsVec = shapePntsRaw(:)';
shapePntsMat = cell2mat(shapePntsVec);
shapePnts = shapePntsMat(~isnan(vecnorm(shapePntsMat,2,2)),:)';
plot3(shapePnts(:,1),shapePnts(:,2),shapePnts(:,3),'sb','MarkerFaceColor','b');axis equal;grid on,grid minor