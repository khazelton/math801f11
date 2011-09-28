function features=ExtractFeatures(basefolder)
% function features=ExtractFeatures(basefolder) extracts several
% features from the images contained in the subfolders of the basefolder.
% basefolder should be '/' terminated.
% The basefolder should contain only subfolders and nothing else.
% The subfolders should contain images in TIF format with extension '.tif'
% They can also contain other data, since the function only processes the
% TIF images.
% The function returns a cell array of all the features that it extracts.
% Each column of the cell features contains the name of the subfolder and a
% cell containing all the relevant data for that subfolder (i.e. the
% extracted midlines, curvatures, average and total curvature, significant
% minima and maxima of the curvature, etc.

%
%  Timing for 10 folders, each with 301 images: 1795 seconds (slow method)
%

MAX_FOLDERS=100; % Process no more than 100 folders
MAX_IMAGES=1000; % Process no more than 1000 images per folder
Cpts=200; % No. of points for curvature computation
cThresh=5; % Threshold for curvatures, everything outside
           % [average/5 , 5*average] is set to 0.
           
Names={'Location','N_images','Midlines','Midsplines','Curvatures',...
       'TotAvgCurvatures','LocalExtrema'};

folders=evalc(['ls ' basefolder]);

subfolders=textscan(folders,'%s',MAX_FOLDERS);

N=size(subfolders{1},1);
if ~N
    error('Empty base folder.');
end

features=cell(1,N);

for i=1:N
    % Process each folder
try
    location=strcat(basefolder, subfolders{1}{i}, '/');
    disp(['Processing '  location]);
    subfolderdata=evalc(['ls ' location]);
    files=textscan(subfolderdata,'%s',MAX_IMAGES);
    if ~size(files{1},1)
        % Empty folder ... go to next one
        continue;
    end
    junk=char(files{1});
    N_img=size(findstr(reshape(junk',1,...
        numel(junk)),'.tif'),2);
    Midlines=cell(1,N_img);
    Midsplines=cell(2,N_img);
    curvatures=cell(1,N_img);
    localextrema=cell(2,N_img);
    minN=100;
    maxN=0;
    for j=1:size(files{1},1)
        pos=findstr(files{1}{j},'.tif');
        if pos
            % Get its number
            imgNum=sscanf(files{1}{j}(1:pos-1),'%d');
            if imgNum < minN
                minN=imgNum;
            end
            if imgNum > maxN
                maxN=imgNum;
            end
        end
    end
    midv=ExtractMidline(location,minN,maxN,...
        [50, .05, -0.5, pi/20, .2, .65, 0 ],1);

    totavgcurvature=zeros(N_img,3);
    for j=1:N_img
        Midlines{j}=midv{j};
        [curve time]=Splines(midv{j});
        Midsplines{1,j}=curve;
        Midsplines{2,j}=time;
        cvs=CurvatureS(curve, time, Cpts, 0);
        curvatures{j}=cvs;
        for k=1:size(cvs,1)-1
            totavgcurvature(j,1)=totavgcurvature(j,1) + ...
                (cvs(k,2)+cvs(k+1,2))*...
                (cvs(k+1,1)-cvs(k,1))/2;
        end
        totavgcurvature(j,3)=cvs(k+1,1);
        totavgcurvature(j,2)=totavgcurvature(j,1)/totavgcurvature(j,3);
        % Adjust the curvatures vector
        cvs(logical(cvs(:,2) < totavgcurvature(j,2)/cThresh),2)=0;
        cvs(logical(cvs(:,2) > totavgcurvature(j,2)*cThresh),2)=0;
        % Now find mins and max
        locex=localminmax(cvs);
        localextrema{1,j}=locex{1};
        localextrema{2,j}=locex{2};       
    end
    % Create the second cell of features
    container=cell(1,7);
    container{1}=location;
    container{2}=N_img;
    container{3}=Midlines;
    container{4}=Midsplines;
    container{5}=curvatures;
    container{6}=totavgcurvature;
    container{7}=localextrema;
    features{i}=cell2struct(container,Names,2);
catch
    container=cell(1,7);
    container{1}=location;
    container{2}=N_img;
    if ~isempty(who('Midlines')) && size(Midlines)
        container{3}=Midlines;
    else
        container{3}=[];
    end
    if ~isempty(who('Midsplines')) && size(Midsplines)
        container{4}=Midsplines;
    else
        container{4}=[];
    end
    if ~isempty(who('curvatures')) && size(curvatures)
        container{5}=curvatures;
    else
        container{5}=[];
    end
    if ~isempty(who('totavgcurvature')) && size(totavgcurvature)
        container{6}=totavgcurvature;
    else
        container{6}=[];
    end
    if ~isempty(who('localextrema')) && size(localextrema)
        container{7}=localextrema;
    else
        container{7}=[];
    end
    features{i}=cell2struct(container,Names,2);
    continue;
end
        
end
