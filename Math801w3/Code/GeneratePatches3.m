function patches=GeneratePatches3(Img, parameters)
% function patches=GeneratePatches3(Img,parameters) generates patches that
% cover the root in the image Img. It does that by first extracting the
% midline, then computing the spline approximation of the midline, finding
% the tangent vectors and then extracting patches with size specified in
% parameters. It rotates patches so that the root is horizontal in each
% patch.
% The parameters are the same as for the midline extraction,
%
%  [ Wsize Delta angThresh dirThresh imThresh endThresh disp_flag]
%
% Patches are spaced according to curvature ...
% TODO: .... Later we might add an extra parameter just for that ...
%
% Default params: [50, .05, -0.5, pi/20, .2, .65, 0 ]
%
% Function returns a structure that contains the number of patches, the
% patches themselves, the times that correspond to the centers of patches,
% the centers of the patches and the angles with respect to the original
% orientation of the patches (the patches are all rotated so that they
% are horizontal).


Names={'N', 'image', 'T', 'center', 'angle', 'Wsize'};

[maxY maxX]=size(Img);

if nargin < 2
    parameters=[50, .05, -0.5, pi/20, .2, .65, 0 ];
end

% Extract the midline
[midc midv]=TraceCurve3(Img, parameters);
% throw away last point
N=size(midv,1);
if N>1
    midv=midv(1:N-1,:);
end

% Take the average window size
Wsize=round(mean(midv(:,3)));

% Make sure it's even
if mod(Wsize,2)
    Wsize=Wsize+1;
end

% Initial patch width
pWidth=Wsize/5;

% Get the spline approximation 
curve=SplinesT(midv(:,1:2));

% Max. number of patches - the max time value is a good approximation to
% the length of the curve, min. patch width is 2
N=round((curve.T(curve.N+1)-curve.T(1))/2);
if N==0
    N=1;
end


content=cell(1,6);

content{2}=cell(N,1);
content{3}=zeros(N,1);
content{4}=zeros(N,2);
content{5}=zeros(N,1);
content{6}=Wsize;

t=curve.T(1);


for i=1:N
    % Are we at the end of the curve?
    if t>curve.T(curve.N+1)
        break;
    end
    [pt vel acc] = EvaluateSplineT(curve, t);

    % Save the point
    content{4}(i,:)=pt;

    pt=round(pt);

    % Save the time and the angle for the i-th patch
    content{3}(i)=t;
    content{5}(i)=atan2(vel(2),vel(1));

    % Now extract the patch, rotate it, crop it and save it
    
    pW=max([Wsize/2, pWidth]);
    if (pt(2) < pW) || (pt(1) < pW)
       % Skip the patch
       continue;
    end

    if (pt(2)+pW > maxX) || (pt(1)+pW > maxY)
       % Skip the patch
       continue;
    end
       
    rowmin=pt(2)-pW + 1;
    rowmax=pt(2)+pW;
    colmin=pt(1)-pW + 1;
    colmax=pt(1)+pW;
    Win=Img(rowmin:rowmax,colmin:colmax);
    WinR=imrotate(Win,content{5}(i)*180/pi,'bicubic');
    pt=size(WinR)/2;
    pW=Wsize*sqrt(2)/3;
    content{2}{i}=WinR(round(pt(1)-pW)+1:round(pt(1)+pW),...
        round(pt(2)-pWidth/2):round(pt(2)+pWidth/2));

    % Estimate the next patch size
    [junk vel2 acc2] = EvaluateSplineT(curve, t+pWidth/norm(vel));
    % Get the average curvature across the patch
    cur=(abs(vel(1)*acc(2)-vel(2)*acc(1))/norm(vel)^3 +...
        abs(vel2(1)*acc2(2)-vel2(2)*acc2(1))/norm(vel2)^3)/2;
    
    % If it's flat, take the largest allowed patch
    if cur==0
        npWidth=Wsize/2;
    else
        npWidth=1/(cur*Wsize);
    end
    
    % Make sure the patch size is between 4 and Wsize
    if npWidth < 4
        npWidth=4;
    end
    if npWidth > Wsize/2
        npWidth=Wsize/2;
    end
    
    % DeltaT - estimate it from the current velocity
    dT=(pWidth+npWidth)/(2*norm(vel)); % (curve.T(curve.N+1)-curve.T(1))/N;
    % Re-estimate the velocity using the next point

    [pt vel2] = EvaluateSplineT(curve, t+dT);
    dT=(pWidth+npWidth)/(norm(vel)+norm(vel2));
    pWidth=npWidth;
    t=t+dT;
end
i=i-1;
if i==0
    i=1;
end
content{1}=i;
content{2}=content{2}(1:i,:);
content{3}=content{3}(1:i,:);
content{4}=content{4}(1:i,:);
content{5}=content{5}(1:i,:);
patches=cell2struct(content, Names, 2);

end







