function fig_hand = draw_body(verts,faces,varargin)
% Inputs:
% verts - [3 x m] matrix, each column is a vertex location.
% faces - [n x 3/4] matrix of faces' vertices. 
%           Entries correspond to a column of verts
% varargin = high_faces -> indices in faces want colored differently
%
% Outputs:
% fig_hand - figure handle

logical_faces_int = false(length(faces),1);

if size(verts,1)~=3
    verts = verts';
end

if nargin == 3
    high_faces = varargin{1};
    logical_faces_int(high_faces) = 1;
end

high_color = [0.87 0.49 0];
body_color = [0.3490 .2392 0];

scrsz = get(0,'ScreenSize');
fig_hand = figure('Position',[1 1 scrsz(3)/2 .6*scrsz(3)]);
hold on; axis equal; set(fig_hand, 'Color', 'w'); axis off
patch('Vertices',verts','Faces',faces(~logical_faces_int,:),'FaceColor',body_color,'Marker','none');
patch('Vertices',verts','Faces',faces(logical_faces_int,:),'FaceColor',high_color,'Marker','none');

end