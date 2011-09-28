function [cubics time]=Splines(midline, disp_flag)
% Function [cubics time]=Splines(midline, disp_flag) returns
% the cubic spline approximation of the curve specified in
% the midline (a NxD matrix), D is the dimension.
% The result is a 1xN-3 cell array of coefficients of 
% cubic polynomials that approximate the curve.
% Each cell contains a 4x2 matrix where the columns correspond
% to x and y coordinate and rows to the coefficients of
% 1, t, t^2 and t^3.
% The corresponding values of the time parameter are returned in
% the variable 'time'
% On interval [time(i), time(i+1)] the approximation in
% cubics{1,i} should be used.
% Also, the last element in time gives a rough estimate of the
% length of the curve.
%
% We approximate the midline using cubic splines and then
% compute the curvature from the splines, using
% k(t)=norm(cross(v(t),a(t)))/norm(v(t))^3
% We use 'arc length' parametrization for knots
% 
% If the midline contains 1, 2 or 3 points, the function returns
% a point, a line, or a quadratic approximation (as a cubic with
% appropriate coefficients set to 0).
% The optional parameter disp_flag determines whether 
% we display the splines as we go along
% default is 0 (no display)


if nargin < 1
    error('Not enough input arguments');
end

if nargin < 2
    disp_flag=0;
end

N=size(midline,1); % Number of points
D=size(midline,2); % Dimension

knots=zeros(1,N+3); % knots are t_2 ... t_N+3, t_2=0
pause_len=0.1;

if N==1
    % Return a point
    time=[0, 0];
    cubics=cell(1,1);
    cubics{1}=zeros(4,D);
    cubics{1}(1,:)=midline(1,:);
    return
elseif N==2
    % Return a line
    time=[0, norm(midline(2,:)-midline(1,:))];
    cubics=cell(1,1);
    cubics{1}=zeros(4,D);
    cubics{1}(1:2,:)=[midline(1,:); ...
        (-midline(1,:)+midline(2,:))/time(2)];  
    return;
elseif N==3
    % Return a quadratic spline
    q1=zeros(3,D);
    cubics=cell(1,1);
    cubics{1}=zeros(4,D);
    time=[0, norm(midline(3,:)-midline(2,:))+...
        norm(midline(2,:)-midline(1,:))];
    l1=[midline(1,:); (-midline(1,:)+midline(2,:))/time(2)];
    l2=[midline(2,:); (-midline(2,:)+midline(3,:))/time(2)];
    q1(1:2,:)=l1;
    q1(2:3,:)=q1(2:3,:) +(-l1+l2)/time(2);
    cubics{1}(1:3,:)=q1;
    return;
end

% N is at least 4 ...
% i=3 case:
knots(1,5)=norm(midline(3,:)-midline(2,:)) ...
            +norm(midline(2,:)-midline(1,:));

for i=4:N-1
    knots(1,i+2)=knots(1,i+1)+norm(midline(i,:)-midline(i-1,:));
end
knots(1,N+1)=knots(1,N+1)+norm(midline(N,:)-midline(N-1,:));
knots(1,N+2)=knots(1,N+1);
knots(1,N+3)=knots(1,N+1);

cubics=cell(1,N-3);

if disp_flag
    N %#ok<NOPRT>
    figure;
    scatter(midline(:,1),midline(:,2),'filled');
    axis equal;
    hold on;
    pause(pause_len);
end

for i=4:N
    % l1, l2, l3 three linear splines
    q1=zeros(3,D); % Two quadratic splines
    q2=zeros(3,D);  % first row for 1, second for t, third for t^2
    c=zeros(4,D); % One cubic spline
    
    l1=[knots(i+1)*midline(i-3,:) - knots(i-2)*midline(i-2,:); ...
        -midline(i-3,:)+midline(i-2,:)];
    l1=l1/(knots(i+1)-knots(i-2));
    l2=[knots(i+2)*midline(i-2,:) - knots(i-1)*midline(i-1,:); ...
        -midline(i-2,:)+midline(i-1,:)];
    l2=l2/(knots(i+2)-knots(i-1));
    l3=[knots(i+3)*midline(i-1,:) - knots(i)*midline(i,:); ...
        -midline(i-1,:)+midline(i,:)];
    l3=l3/(knots(i+3)-knots(i));
    if disp_flag && D==2
        t=[knots(i-2);knots(i+1)];
        data=[ones(size(t,1),1),t]*l1;
        plot(data(:,1),data(:,2),'r');
        t=[knots(i-1);knots(i+2)];
        data=[ones(size(t,1),1),t]*l2;
        plot(data(:,1),data(:,2),'r');
        t=[knots(i);knots(i+3)];
        data=[ones(size(t,1),1),t]*l3;
        plot(data(:,1),data(:,2),'r');
        pause(pause_len);
    end

    q1(1:2,:)=knots(i+1)*l1 - knots(i-1)*l2;
    q1(2:3,:)=q1(2:3,:) -l1+l2;    
    q1=q1/(knots(i+1)-knots(i-1));
    
    q2(1:2,:)=knots(i+2)*l2 - knots(i)*l3;
    q2(2:3,:)=q2(2:3,:) -l2+l3;    
    q2=q2/(knots(i+2)-knots(i));

    if disp_flag && D==2
        dt=(knots(i+1)-knots(i-1))/50;
        t=(knots(i-1):dt:knots(i+1))';
        data=[ones(size(t,1),1),t,t.^2]*q1;
        plot(data(:,1),data(:,2),'g');
        dt=(knots(i+2)-knots(i))/50;
        t=(knots(i):dt:knots(i+2))';
        data=[ones(size(t,1),1),t,t.^2]*q2;
        plot(data(:,1),data(:,2),'g');
        pause(pause_len);
    end
    
    c(1:3,:)=knots(i+1)*q1 - knots(i)*q2;
    c(2:4,:)=c(2:4,:) - q1+q2;
    c=c/(knots(i+1)-knots(i));
    
    if disp_flag && D==2
        dt=(knots(i+1)-knots(i))/50;
        t=(knots(i):dt:knots(i+1))';
        data=[ones(size(t,1),1),t,t.^2,t.^3]*c;
        plot(data(:,1),data(:,2),'b','linewidth',2);
        scatter([1 knots(i) knots(i)^2 knots(i)^3]*c(:,1), ...
            [1 knots(i) knots(i)^2 knots(i)^3]*c(:,2),'r','filled');
        scatter([1 knots(i+1) knots(i+1)^2 knots(i+1)^3]*c(:,1), ...
            [1 knots(i+1) knots(i+1)^2 knots(i+1)^3]*c(:,2),'g','filled');
        pause(pause_len);
    end
    
    cubics{1,i-3}=c;
end
time=knots(4:N+1);

end