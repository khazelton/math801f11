function curvatures=CurvatureS(curve, time, N, disp_flag)
% Function curvatures=CurvatureS(curve, time, N) returns the 
% curvature along the curve specified by the cell array of
% cubic spline approximations with values of time parameter in time
% at N points along the curve. It returns triplets
% (arclength, curvature, t) along the curve.
% The arclenght represents the length of the curve up to that point,
% curvature is the curvature at that point and t is the time at which
% this point appears on the curve. This allows for easier visualization.
% It assumes that the interpolation between time(i)
% and time(i+1) is given in curve{1,i}.
% time should be a row vector with exactly one element more
% than the curve cell array
%
% The curvature is computed using
% k(t)=norm(cross(v(t),a(t)))/norm(v(t))^3


if nargin < 3
    error('Not enough input arguments');
end

if nargin < 4
    disp_flag=0;
end

M=size(curve,2);
if size(time,2) ~= M+1
    error('Data size mismatch');
end

maxT=time(M+1);


dt=(maxT-time(1))/N;

curvatures=zeros(N,3);

% Internal parameters ...
NN = 10; % No. of in-between points

t=time(1); % initial time
curcell=1; % current cell

if disp_flag
    figure;
    hold on;
    axis equal;
end;

len=0; %Initial position

for i=1:N
    while ( curcell <= M) && (t > time(curcell+1))
        curcell=curcell+1;
    end
    if curcell==M+1
        break;
    end
    ptpoly=[1 t t^2 t^3];
    velpoly=[0 1 2*t 3*t^2];
    accpoly=[0 0 2 6*t];
    vel=velpoly*curve{1,curcell};
    acc=accpoly*curve{1,curcell};
    if disp_flag
        pt=ptpoly*curve{1,curcell};
        scatter(pt(1),pt(2),'b','filled');
        line([pt(1); pt(1)+vel(1)],[pt(2); pt(2)+vel(2)],'Color','g');
        line([pt(1); pt(1)+acc(1)],[pt(2); pt(2)+acc(2)],'Color','r');
        pause(0.1);
    end
    curvatures(i,:)=[len, norm(cross([vel, 0],[acc,0]))/norm(vel)^3, t];
    % Now compute the new length
    dtt=dt/NN;
    tt=t+dtt;
    prevcell=curcell;
    for j=1:NN
        while (curcell <= M) && (tt > time(curcell+1))
            curcell=curcell+1;
        end
        if curcell==M+1
            break;
        end
        nptpoly=[1 tt tt^2 tt^3];
        len=len+norm(nptpoly*curve{1,curcell} - ptpoly*curve{1,prevcell});
        ptpoly=nptpoly;
        tt=tt+dtt;
        prevcell=curcell;
    end
    t=t+dt;
end



end