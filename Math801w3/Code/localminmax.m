function loc_extrema=localminmax(data)
% function loc_extrema=localminmax(data) finds local extrema from the
% matrix data (a Nx2 matrix, first column contains positions, second
% contains the function values. The function returns a cell with two
% elements, the first contains local maxima, the second contains local
% minima, both as vectors with positions of extrema within data.
% data should have all distinct entries in the first column.
% If data has more than 2 rows we return only the local extrema inside the
% interval defined by the first column and ignore the endpoints.

% Ignore columns 3 and higher
data=data(:,1:2);

N=size(data,1);
loc_extrema=cell(1,2);

if N==1
    loc_extrema{1}=1;
    loc_extrema{2}=1;
    return
elseif N==2
    if data(1,2) >= data(2,2)
        loc_extrema{1}=1;
        loc_extrema{2}=2;
    else
        loc_extrema{1}=2;
        loc_extrema{2}=1;
    end
    return
end

loc_extrema{1}=zeros(N-2,1);
loc_extrema{2}=zeros(N-2,1);

minpos=0;
maxpos=0;

dxy1=data(2,:)-data(1,:);

if dxy1(1) <= 0
    disp('Bad data - non-increasing x entries.')
    loc_extrema{1}=[];
    loc_extrema{2}=[];
    return
end

df1=dxy1(2)/dxy1(1);

for i=2:N-1
    dxy2=data(i+1,:)-data(i,:);
    if dxy2(1) <= 0
        disp('Bad data - non-increasing x entries.')
        break;
    end
    df2=dxy2(2)/dxy2(1);
    if df1*df2 < 0
        % Derivative changes sign at i-th entry ...
        if df1 < 0
            % \/ -> min
            minpos=minpos+1;
            loc_extrema{2}(minpos,:)=i;
        else
            % /\ -> max
            maxpos=maxpos+1;
            loc_extrema{1}(maxpos,:)=i;
        end
    end
    df1=df2;
end

if maxpos
    loc_extrema{1}=loc_extrema{1}(1:maxpos);
else
    loc_extrema{1}=[];
end
if minpos
    loc_extrema{2}=loc_extrema{2}(1:minpos);
else
        loc_extrema{2}=[];
end

end
        

    
    
    
    
    
    
