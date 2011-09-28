function [Coords vertices]=TraceCurve3(Img, parameters, slow)
% Do SVD for windows going over the image
% [Coords vertices]=TraceCurve(Img, init_dir) returns
% coordinates of windows that trace over the curve given
% in Img. In addition, it returns the vertices of the graph that
% fits the image together with the width of the curve at each vertex.
% The initial direction specifies whether we 
% search horizontally or vertically, init_dir=0 for vertical,
% 1 for horizontal.
%
% The parameters are
% 
% [ Wsize Delta angThresh dirThresh imThresh endThresh disp_flag]
%
% where Wsize is the initial windows size (default: 1/8 of the smaller
% dimension of the image);
%
% Delta is the fraction of Wsize we move each time, default 1/10;
%
% angThreshold is the angle threshold, we don't allow to go in directions
% where the cos of the angle with the previous one is <angThresh, default
% is -.1, i.e. direction changes of more than 95 degrees are not allowed
%
% dirThresh (in radians): we measure cummulative angle change as we move 
% along the curve, if it exceeds the dirThresh we issue a new vertex.
% The default value is pi/15;
%
% imThresh is the pixel cutoff as a fraction of maximum, we ignore all
% pixels that are below imThresh*max;
%
% endThresh decides when to stop; at the end of the tube the singular
% values get close to each other; default is .6 (as proportion of
% the largest sing. value to the sum of all
%
% disp_flag determines whether the progress is graphically displayed,
% 1 means yes, 0 (default) means no.
%
% The last optional flag 'slow' indicates whether we do fast tracing
% (slow=0) or slower and more accurate tracing (slow=1) where we
% recenter data by masking. Default is slow=0.

if nargin < 1
    error('Not enough input arguments')
end

[maxy maxx] = size(Img);

if nargin < 2
    parameters=[ round(min(maxx,maxy)/8), .1, -.1, pi/15, .2, .6, 0 ];
end

if nargin < 3
    slow=0;
end

if length(parameters) < 7
    error('Not enough parameters');
end

Wsize=parameters(1);
if mod(Wsize,2)
    % Make it even
    Wsize=Wsize+1;
end
Delta=max(Wsize*parameters(2),2);
angThreshold=parameters(3);
dirThreshold=parameters(4);
endThresh=parameters(6);
disp_flag=parameters(7);

data=double(max(max(Img)))-double(Img);

if slow
    % replicate the first column Wsize/2 times so that we get the
    % trace from the beginning of the curve
    iWsize=Wsize/2; % Record is so that we can subtract it at the end
    data=[data(:,1)*ones(1,iWsize),data];
end

maxData=max(max(data));
imThreshold=parameters(5)*maxData;
maxVert=max(maxx,maxy);


Coords=zeros(round(((maxx+maxy)/(2*Delta))^2),3);

tubeWFactor=2; % Default: 2


best_i=1;
best_j=1;
best_norm=0;
% Do vertical search of initial window
for i=1:round(Delta):(maxy-Wsize)
    Win=data(i:i+Wsize-1,1:Wsize);
    cur_norm=norm(Win,'fro');
    if cur_norm>best_norm
        best_norm=cur_norm;
        best_i=i;
    end;
end;

slice=1; % no. of slices (steps)
pause_len=.01; % Length of pause for display
cangrow=1; % Can the tube grow?
cangrowFactor=4; % How many times can it grow between forks?

Win=data(best_i:best_i+Wsize-1,best_j:best_j+Wsize-1);
% Recenter data and adjust window size and Delta
[ppp2 ppp1]=find(Win>imThreshold);
N=length(ppp1);
pairs=[ppp1, ppp2];
MEAN=mean(pairs,1);
pairs_N=pairs-(MEAN(ones(N,1),:));
[pU pS]=svd(pairs_N'*pairs_N/N);
oldWsize=Wsize;
Wsize=round(tubeWFactor*(-1+sqrt(1+12*pS(2,2))));
if mod(Wsize,2)
    Wsize=Wsize+1;
end;
Coords(slice,3)=Wsize;


if disp_flag
    if oldWsize ~= Wsize
        Wsize %#ok<NOPRT>
    end
end

best_i=best_i+round(MEAN(1,2));
best_j=Wsize/2;

Win=data(best_i-Wsize/2+1:best_i+Wsize/2, ...
             best_j-Wsize/2+1:best_j+Wsize/2);

if disp_flag
    figure('Position',[64, 10, min(maxx+100,1500),...
                                    min(maxy+100,1000)]);
    imshow(Img,'InitialMagnification','fit');
    hold on
    if slow
        h=squ(best_j, best_i, Wsize);
        scatter(best_j-iWsize,best_i,'g')
    else
        scatter(best_j,best_i,'g')
    end
    pause(pause_len);
end

vertices=zeros(maxVert,3);
vertices(1,:)=[0, 0, 0];
vertices(2,:)=[best_j, best_i, Wsize ];
Nvert=2;

old_dir=[1;0];
oldbest_i=0;
oldbest_j=0;
direct=1;
rep=0;
cum_dir_ch=0;

while best_j<maxx-Wsize/2
    % Form the pairs that work
    [ppp2 ppp1]=find(Win>imThreshold);
    N=length(ppp1);
    pairs=[ppp1, ppp2];

    % normalize
    MEAN=mean(pairs,1);
    pairs_N=pairs-(MEAN(ones(N,1),:));
    [pU pS]=svd(pairs_N'*pairs_N/N);
   
    angle_ch=old_dir'*pU(:,1);
    if angle_ch < angThreshold
        % bad dir, one is allowed, but if it keeps happening bail out
        if rep<5
            direct=-direct;
            angle_ch=-angle_ch;
            rep=rep+1;
        else
            disp('Bad direction.');
            break;
        end;
    else
        rep=0;
    end;
    % Do we need to issue a new vertex?
    angle_ch=acos(angle_ch);
    cum_dir_ch = cum_dir_ch + angle_ch;
    if cum_dir_ch > dirThreshold
        % Issue new vertex
        cum_dir_ch=0;
        Nvert=Nvert+1;
        vertices(Nvert,:)=[best_j, best_i, Wsize];
        if disp_flag
            if slow
                scatter(best_j-iWsize, best_i,'r','filled');
            else
                scatter(best_j, best_i,'r','filled');
            end;
            pause(pause_len);
        end
    end;

    Coords(slice,1:2)=[best_j, best_i];
    % Assign new coordinates and recenter data ...
    
    tubediam=(-1+sqrt(1+12*pS(2,2)));
    
    
    if slow
        MEAN2=[0,0];
        k=0;
        vecU=direct*pU(:,1);
        lowerThr = MEAN*vecU-tubediam/3;
        %lowerThr = [Wsize/2, Wsize/2]*vecU-tubediam/3;
        upperThr = lowerThr+2*tubediam/3;

        if vecU(1) > 0
            for i=1:Wsize
                proj=i*vecU(2);
                for j=1:Wsize
                    if (Win(i,j) > imThreshold)
                        proj2=proj+j*vecU(1);
                        if proj2 >= upperThr
                            break;
                        end
                        if proj2>lowerThr
                            k=k+1;       
                            MEAN2=MEAN2+[j,i];
                        end
                    end
                end
            end
        else
            for i=1:Wsize
                proj=i*vecU(2);
                for j=1:Wsize
                    if (Win(i,j) > imThreshold)
                        proj2=proj+j*vecU(1);
                        if proj2 <= lowerThr
                            break;
                        end
                        if proj2<upperThr
                            k=k+1;       
                            MEAN2=MEAN2+[j,i];
                        end
                    end
                end
            end
        end
        if k==0
            MEAN2=[Wsize/2, Wsize/2];
        else
            MEAN2=MEAN2/k;
        end
    else
        MEAN2=MEAN;
    end
    
    correction=([MEAN2(1,1)-Wsize/2,MEAN2(1,2)-Wsize/2]*pU(:,2))*pU(:,2);
    best_j=best_j + round(direct*pU(1,1)*Delta + correction(1));
    best_i=best_i + round(direct*pU(2,1)*Delta + correction(2));

    if pS(1,1)/trace(pS) < endThresh
        % at the end ...
        Nvert=Nvert+1;
        vertices(Nvert,:)=[best_j,best_i, Wsize];
        if disp_flag
            if slow
                scatter(best_j-iWsize, best_i,'r','filled');
            else
                scatter(best_j, best_i,'r','filled');
            end
            pause(pause_len);
        end;
        break;
    end;

    if ((oldbest_i==best_i) && (oldbest_j==best_j))
        %same point, issue a vertex and stop
        % This should never happen.
        Nvert=Nvert+1;
        vertices(Nvert,:)=[best_j,best_i, Wsize];
        if disp_flag
            if slow
                scatter(best_j-iWsize, best_i,'r','filled');
            else
                scatter(best_j, best_i,'r','filled');
            end
            pause(pause_len);
        end;
        break;
    end
    oldbest_i=best_i;
    oldbest_j=best_j;

    % Adjust window size
    oldWsize=Wsize;
    Wsize=round(tubeWFactor*tubediam);
    if mod(Wsize,2)
        Wsize=Wsize+1;
    end
    % Do we grow too much?
    if Wsize>oldWsize
        if cangrow
            Wsize=oldWsize+2;
            % If we grew more than allowed, don't grow
            % until the next fork
            cangrow=cangrow-1;
        else                
            Wsize=oldWsize;
        end
    elseif Wsize < oldWsize
        % We shrink, so allow growth later on
        cangrow=min(cangrow+1,cangrowFactor);
    end
    Coords(slice,3)=Wsize;
    Delta=max(Wsize*parameters(2),2);
    if disp_flag
        if oldWsize ~= Wsize
            Wsize %#ok<NOPRT>
        end
    end
    if best_j<=Wsize/2
        best_j=Wsize/2+1;
    end
    nextWin=data(best_i-Wsize/2:best_i+Wsize/2-1,...
        best_j-Wsize/2:best_j+Wsize/2-1);

    slice=slice+1;
    if disp_flag
        if slow
            scatter(best_j-iWsize, best_i, 'ro');
        else
            scatter(best_j, best_i, 'ro');
        end
        slice %#ok<NOPRT>
        pause(pause_len);
    end
    old_dir=pU(:,1);
    Win=nextWin;
    if slice> round((maxx/Delta)^2)
        disp('Too many slices');
        break;
    end;

end;
    
if Nvert>2
    if vertices(2,:) == vertices(3,:)
        vertices=vertices(3:Nvert,:);
    else
        vertices=vertices(2:Nvert,:);
    end    
end
Coords(slice,:)=[best_j, best_i, Wsize];
Coords=Coords(1:slice,:);
if slow
    Coords(:,1)=Coords(:,1)-iWsize;
    vertices(:,1)=vertices(:,1)-iWsize;
end

end

function h=squ(x, y, len)
   h=line([x-len/2,x+len/2-1,x+len/2-1,x-len/2,x-len/2],...
       [y-len/2,y-len/2,y+len/2-1,y+len/2-1,y-len/2],'color','b');
end
