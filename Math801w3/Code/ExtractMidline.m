function [midverts midlines]=ExtractMidline(folder, n_b, n_e, parameters, slow)
% Function 
% [midverts midlines]=ExtractMidline(folder, n_b, n_e, parameters, slow) 
% extracts midlines from the tif images found in the folder 'folder'. The 
% images are assumed to be named n_b.tif, n_b+1.tif, ... n_e.tif,
% where n_b and n_e are integers, usually 0 and 300.
% The extracted midlines are placed in the cell midlines.
% The extracted vertices are placed in the cell midverts.
% Parameters have the same meaning as for TraceCurve3:
%  
% [ Wsize Delta angThresh dirThresh imThresh endThresh disp_flag]
%
% If slow=1, the precise (and slower) midline tracing is used in
% TraceCurve3, otherwise the fast one is used (default)

if nargin < 4
    parameters=[50, .05, -0.5, pi/15, .2, .65, 0 ];
end;

if nargin < 5
    slow=0;
end

midlines=cell(1,n_e-n_b+1);
midverts=cell(1,n_e-n_b+1);

for i=n_b:n_e
    I=imread(strcat(folder,int2str(i),'.tif'));
    [coords verts]=TraceCurve3(I,parameters,slow);
    midlines{1,i-n_b+1}=coords(:,1:2);
    midverts{1,i-n_b+1}=verts(:,1:2);
    if parameters(7) || slow
        % display flag is set or we are doing slow processing ...
        if mod(i,10)==0
            i %#ok<NOPRT>
        end
    end
end

end
    

    