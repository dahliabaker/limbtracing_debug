clear;clc;
%%
% settings
computeConvexHull = 0;

% set path to image directory
% imDir = '.\0bennu_sim_IR\bennu_automated_images\';% bennu
imDir = '.\0itokawa_sim_IR\itokawa_automated_images\';% itokawa

% compute number of images in directory
numFiles = length(dir([imDir,'render*']));

% set figure size
h = figure(10);
h.Position = [70,0,1000,500];

% initialize arrays
irSilhoutte = cell(numFiles,2);


% start processing images
for idxIm = 1:numFiles
    %% finding silhoutte
    % read image
    im = imread([imDir,'render',num2str(idxIm),'.png']);
    
    % convert to grayscale
    im = rgb2gray(im);
    
    % remove random background noise
    im(im<uint8(2)) = uint8(0);
    
    % compute image silhouette
    imthresh = im*1000;
    
    imSilhoutte = zeros(size(im),'uint8');
    for ii = 1:size(im,1)
        if any(im(ii,:)>0)
            c = find(im(ii,:)>0,1,'first');
            imSilhoutte(ii,c) = uint8(1000);
            c = find(im(ii,:)>0,1,'last');
            imSilhoutte(ii,c) = uint8(1000);
        end
    end
    for ii = 1:size(im,2)
        if any(im(:,ii)>0)
            c = find(im(:,ii)>0,1,'first');
            imSilhoutte(c,ii) = uint8(1000);
            c = find(im(:,ii)>0,1,'last');
            imSilhoutte(c,ii) = uint8(1000);
        end
    end
    
    % set figure title
    sgtitle(['render',num2str(idxIm)])
    
    % find silhoutte
    [c,r] = find(imSilhoutte);
    % optional convex hull computation
    if computeConvexHull
        DT = delaunayTriangulation(r,c);
        
        C = convexHull(DT);
    end
    
    %% plotting
    ax1 = subplot(131);
    cla
    
    % plot original image
    imshow(im)
    colormap(gray);
    axis square

    % plot image gradient
    ax2 = subplot(132);
    cla
    imshow(imSilhoutte)
    colormap(gray)
    axis square
    
    % plot silhoutte (and optionally the convex hull)
    ax3 = subplot(133);
    cla
    plot(r,c,'.')
    
    % optionally plot convex hull
    if computeConvexHull
        hold on
        plot(DT.Points(C,2),DT.Points(C,1),'r+-','MarkerSize',10)
    end
    axis square
    axis([0,1024,0,1024])
    set(ax3,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    drawnow
    
    % save data
    irSilhoutte{idxIm,1} = c;
    irSilhoutte{idxIm,2} = r;
    
end