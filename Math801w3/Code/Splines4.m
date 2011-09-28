function [quartics time]=Splines4(midline, disp_flag)
% Function [quartics time]=Splines4(midline, disp_flag) returns
% the quartic spline approximation of the curve specified in
% the midline (a NxD matrix), D is the dimension.
% The result is a 1xN-3 cell array of coefficients of 
% quartic polynomials that approximate the curve.
% Each cell contains a 4x2 matrix where the columns correspond
% to x and y coordinate and rows to the coefficients of
% 1, t, t^2 and t^3.
% The corresponding values of the time parameter are returned in
% the variable 'time'
% On interval [time(i), time(i+1)] the approximation in
% quartics{1,i} should be used.
% Also, the last element in time gives a rough estimate of the
% length of the curve.
%
% We approximate the midline using quartic splines and then
% We use 'arc length' parametrization for knots
% 
% If the midline contains 1, 2, 3 or 4 points, the function returns
% a point, a line, a quadratic approximation or a cubic approximation 
% (as a quartic with appropriate coefficients set to 0).
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

knots=zeros(1,N+4); % knots are t_2 ... t_N+3, t_2=0
pause_len=0.1;

if disp_flag && D==2
    N %#ok<NOPRT>
    figure;
    scatter(midline(:,1),midline(:,2),'filled');
    axis equal;
    hold on;
    pause(pause_len);
end

if N==1
    % Return a point
    time=[0, 0];
    quartics=cell(1,1);
    quartics{1}=zeros(5,D);
    quartics{1}(1,:)=midline(1,:);
    return
elseif N==2
    % Return a line
    time=[0, norm(midline(2,:)-midline(1,:))];
    quartics=cell(1,1);
    quartics{1}=zeros(5,D);
    quartics{1}(1:2,:)=[midline(1,:); ...
        (-midline(1,:)+midline(2,:))/time(2)];
    if disp_flag && D==2
        data=[ones(size(time,2),1),time']*quartics{1}(1:2,:);
        plot(data(:,1),data(:,2),'r');
        pause(pause_len);
    end
    return;
elseif N==3
    % Return a quadratic spline
    q1=zeros(3,D);
    quartics=cell(1,1);
    quartics{1}=zeros(5,D);
    time=[0, norm(midline(3,:)-midline(2,:))+...
        norm(midline(2,:)-midline(1,:))];
    l1=[midline(1,:); (-midline(1,:)+midline(2,:))/time(2)];
    l2=[midline(2,:); (-midline(2,:)+midline(3,:))/time(2)];
    q1(1:2,:)=l1;
    q1(2:3,:)=q1(2:3,:) +(-l1+l2)/time(2);
    quartics{1}(1:3,:)=q1;
    if disp_flag && D==2
        data=[ones(size(time,2),1),time']*l1;
        plot(data(:,1),data(:,2),'r');
        data=[ones(size(time,2),1),time']*l2;
        plot(data(:,1),data(:,2),'r');
        dt=time(2)/50;
        t=(0:dt:time(2))';
        data=[ones(size(t,1),1),t,t.^2]*q1;
        plot(data(:,1),data(:,2),'g');
        pause(pause_len);
    end
    return;
elseif N==4
    % Return a cubic spline
    quartics=cell(1,1);
    quartics{1}=zeros(5,D);
    
    time=zeros(1,2);
    time(2)=norm(midline(2,:)-midline(1,:))+...
        norm(midline(3,:)-midline(2,:))+...
        norm(midline(4,:)-midline(3,:));
    
    q1=zeros(3,D); % Two quadratic splines
    q2=zeros(3,D);  % first row for 1, second for t, third for t^2
    c=zeros(4,D); % One cubic spline
    
    l1=[midline(1,:); (-midline(1,:)+midline(2,:))/time(2)];
    l2=[midline(2,:); (-midline(2,:)+midline(3,:))/time(2)];
    l3=[midline(3,:); (-midline(3,:)+midline(4,:))/time(2)];

    q1(1:2,:)=l1;
    q1(2:3,:)=q1(2:3,:)+(l2-l1)/time(2);
    
    q2(1:2,:)=l2;
    q2(2:3,:)=q2(2:3,:) +(l3-l2)/time(2);

    c(1:3,:)=q1;
    c(2:4,:)=c(2:4,:) +(q2-q1)/time(2);
    
    quartics{1}(1:4,:)=c;
    
    if disp_flag && D==2
        data=[ones(size(time,2),1),time']*l1;
        plot(data(:,1),data(:,2),'r');
        data=[ones(size(time,2),1),time']*l2;
        plot(data(:,1),data(:,2),'r');
        data=[ones(size(time,2),1),time']*l3;
        plot(data(:,1),data(:,2),'r');
        dt=time(2)/50;
        t=(0:dt:time(2))';
        data=[ones(size(t,1),1),t,t.^2]*q1;
        plot(data(:,1),data(:,2),'g');
        data=[ones(size(t,1),1),t,t.^2]*q2;
        plot(data(:,1),data(:,2),'g');
        data=[ones(size(t,1),1),t,t.^2,t.^3]*c;
        plot(data(:,1),data(:,2),'b');
        pause(pause_len);
    end
    return;
end

% N is at least 5 ...
% i=4 case:
knots(1,6)=norm(midline(4,:)-midline(3,:))+...
            norm(midline(3,:)-midline(2,:))+...
            norm(midline(2,:)-midline(1,:));
        
for i=5:N-1
    knots(1,i+2)=knots(1,i+1)+norm(midline(i,:)-midline(i-1,:));
end
knots(1,N+1)=knots(1,N+1)+norm(midline(N,:)-midline(N-1,:));
knots(1,N+2)=knots(1,N+1);
knots(1,N+3)=knots(1,N+1);
knots(1,N+4)=knots(1,N+1);

quartics=cell(1,N-4);



for i=5:N
    % l1, l2, l3, l4 four linear splines
    q1=zeros(3,D); % Three quadratic splines
    q2=zeros(3,D);  % first row for 1, second for t, third for t^2
    q3=zeros(3,D);
    c1=zeros(4,D); % Two cubic splines
    c2=zeros(4,D);
    k=zeros(5,D); % One quartic spline
    
    
    l1=[knots(i+1)*midline(i-4,:) - knots(i-3)*midline(i-3,:); ...
        -midline(i-4,:)+midline(i-3,:)];
    l1=l1/(knots(i+1)-knots(i-3));
    l2=[knots(i+2)*midline(i-3,:) - knots(i-2)*midline(i-2,:); ...
        -midline(i-3,:)+midline(i-2,:)];
    l2=l2/(knots(i+2)-knots(i-2));
    l3=[knots(i+3)*midline(i-2,:) - knots(i-1)*midline(i-1,:); ...
        -midline(i-2,:)+midline(i-1,:)];
    l3=l3/(knots(i+3)-knots(i-1));
    l4=[knots(i+4)*midline(i-1,:) - knots(i)*midline(i,:); ...
        -midline(i-1,:)+midline(i,:)];
    l4=l4/(knots(i+4)-knots(i));
    if disp_flag && D==2
        t=[knots(i-3);knots(i+1)];
        data=[ones(size(t,1),1),t]*l1;
        plot(data(:,1),data(:,2),'r');
        t=[knots(i-2);knots(i+2)];
        data=[ones(size(t,1),1),t]*l2;
        plot(data(:,1),data(:,2),'r');
        t=[knots(i-1);knots(i+3)];
        data=[ones(size(t,1),1),t]*l3;
        plot(data(:,1),data(:,2),'r');
        t=[knots(i);knots(i+4)];
        data=[ones(size(t,1),1),t]*l4;
        plot(data(:,1),data(:,2),'r');
        pause(pause_len);
    end

    q1(1:2,:)=knots(i+1)*l1 - knots(i-2)*l2;
    q1(2:3,:)=q1(2:3,:) -l1+l2;    
    q1=q1/(knots(i+1)-knots(i-2));
    
    q2(1:2,:)=knots(i+2)*l2 - knots(i-1)*l3;
    q2(2:3,:)=q2(2:3,:) -l2+l3;    
    q2=q2/(knots(i+2)-knots(i-1));
    
    q3(1:2,:)=knots(i+3)*l3 - knots(i)*l4;
    q3(2:3,:)=q3(2:3,:) -l3+l4;    
    q3=q3/(knots(i+3)-knots(i));

    if disp_flag && D==2
        dt=(knots(i+1)-knots(i-2))/50;
        t=(knots(i-2):dt:knots(i+1))';
        data=[ones(size(t,1),1),t,t.^2]*q1;
        plot(data(:,1),data(:,2),'g');
        dt=(knots(i+2)-knots(i-1))/50;
        t=(knots(i-1):dt:knots(i+2))';
        data=[ones(size(t,1),1),t,t.^2]*q2;
        plot(data(:,1),data(:,2),'g');
        dt=(knots(i+3)-knots(i))/50;
        t=(knots(i):dt:knots(i+3))';
        data=[ones(size(t,1),1),t,t.^2]*q3;
        plot(data(:,1),data(:,2),'g');
        pause(pause_len);
    end
    
    c1(1:3,:)=knots(i+1)*q1 - knots(i-1)*q2;
    c1(2:4,:)=c1(2:4,:) - q1+q2;
    c1=c1/(knots(i+1)-knots(i-1));
    c2(1:3,:)=knots(i+2)*q2 - knots(i)*q3;
    c2(2:4,:)=c2(2:4,:) - q2+q3;
    c2=c2/(knots(i+2)-knots(i));
    
    if disp_flag && D==2
        dt=(knots(i+1)-knots(i-1))/50;
        t=(knots(i-1):dt:knots(i+1))';
        data=[ones(size(t,1),1),t,t.^2,t.^3]*c1;
        plot(data(:,1),data(:,2),'b');
        dt=(knots(i+2)-knots(i))/50;
        t=(knots(i):dt:knots(i+2))';
        data=[ones(size(t,1),1),t,t.^2,t.^3]*c2;
        plot(data(:,1),data(:,2),'b'); 
    end
    
    k(1:4,:)=knots(i+1)*c1 - knots(i)*c2;
    k(2:5,:)=k(2:5,:) - c1 + c2;
    k=k/(knots(i+1)-knots(i));
    
    if disp_flag && D==2
        dt=(knots(i+1)-knots(i))/50;
        t=(knots(i):dt:knots(i+1))';
        data=[ones(size(t,1),1),t,t.^2,t.^3,t.^4]*k;
        plot(data(:,1),data(:,2),'color',[1,.5,0],'linewidth',2);
        scatter([1 knots(i) knots(i)^2 knots(i)^3 knots(i)^4]*k(:,1), ...
            [1 knots(i) knots(i)^2 knots(i)^3 knots(i)^4]*k(:,2),'r','filled');
        scatter([1 knots(i+1) knots(i+1)^2 knots(i+1)^3 knots(i+1)^4]*k(:,1), ...
            [1 knots(i+1) knots(i+1)^2 knots(i+1)^3 knots(i+1)^4]*k(:,2),'g','filled');
        pause(pause_len);
    end
    
    quartics{1,i-4}=k;
end
time=knots(5:N+1);

end