function patches=GeneratePatches(Img, parameters)
% function patches=GeneratePatches(Img) generates patches that cover
% the root in the image Img. It does that by first extracting the midline,
% then computing the spline approximation of the midline, finding the
% tangent vectors and then extracting patches with size specified in
% parameters. It rotates patches so that the root is horizontal in each
% patch.
% The parameters are the same as for the midline extraction,
%
%  [ Wsize Delta angThresh dirThresh imThresh endThresh disp_flag]
%
% Patches are spaced by Wsize/5 
% TODO: .... Later we might add an extra parameter just for that ...
%
% Default params: [50, .05, -0.5, pi/20, .2, .65, 0 ]
%
% Function returns a structure that contains the number of patches, the
% patches themselves, the times that correspond to the centers of patches,
% the centers of the patches and the angles with respect to the original
% orientation of the patches (the patches are all rotated so that they
% are horizontal).


Names={'N', 'image', 'T', 'center', 'angle'};

[maxY maxX]=size(Img);

if nargin < 2
    parameters=[50, .05, -0.5, pi/20, .2, .65, 0 ];
end

% Extract the midline
[midc midv]=TraceCurve3(Img, parameters);

% Take the average window size
Wsize=round(mean(midv(:,3)));

% Make sure it's even
if mod(Wsize,2)
    Wsize=Wsize+1;
end


% Get the spline approximation
curve=SplinesT(midv(:,1:2));

% Number of patches - the max time value is a good approximation to the
% length of the curve
N=round((curve.T(curve.N+1)-curve.T(1))*5/Wsize);
if N==0
    N=1;
end


content=cell(1,5);
content{1}=N;
content{2}=cell(N,1);
content{3}=zeros(N,1);
content{4}=zeros(N,2);
content{5}=zeros(N,1);


t=curve.T(1);


for i=1:N
   [pt vel] = EvaluateSplineT(curve, t);
   
   % Save the point
   content{4}(i,:)=pt;
   
   pt=round(pt);
   
   % Save the time and the angle for the i-th patch
   content{3}(i)=t;
   content{5}(i)=atan2(vel(2),vel(1));
   
   % Now extract the patch, rotate it, crop it and save it
   
   if (pt(2) < Wsize/2) || (pt(1) < Wsize/2)
       % Skip the patch
       continue;
   end
   
   if (pt(2)+Wsize/2 > maxX) || (pt(1)+Wsize/2 > maxY)
       % Skip the patch
       continue;
   end
   
   rowmin=pt(2)-Wsize/2 + 1;
   rowmax=pt(2)+Wsize/2;
   colmin=pt(1)-Wsize/2 + 1;
   colmax=pt(1)+Wsize/2;
   Win=Img(rowmin:rowmax,colmin:colmax);
   WinR=imrotate(Win,content{5}(i)*180/pi,'bicubic');
   pt=size(WinR)/2;
   % sqrt(24)/10 = 0.4899
   content{2}{i}=WinR(round(pt(1)-Wsize*0.49)+1:round(pt(1)+Wsize*0.49), ...
       round(pt(2)-Wsize/10):round(pt(2)+Wsize/10));
   
   % DeltaT - estimate it from the current velocity
   dT=Wsize/(5*norm(vel)); % (curve.T(curve.N+1)-curve.T(1))/N;
   t=t+dT;
end

patches=cell2struct(content, Names, 2);

end





