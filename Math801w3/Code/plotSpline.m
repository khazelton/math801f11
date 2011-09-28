function [points len]=plotSpline(curve, time, N)
% Function [points len]=plotSpline(curve, times, N) computes
% N points on the cubic spline interpolation specified
% in the cell curve, corresponding to the time values.
% It assumes that the interpolation between time(i)
% and time(i+1) is given in curve{1,i}.
% time should be a row vector with exactly one element more
% than the curve cell array

if nargin < 3
    error('Not enough input arguments');
end

M=size(curve,2);
if size(time,2) ~= M+1
    error('Data size mismatch');
end

maxT=time(M+1);


dt=(maxT-time(1))/N;

points=zeros(N,2);

t=time(1); % initial time
curcell=1; % current cell

len=0;

for i=1:N
    while (t > time(curcell+1)) && ( curcell <= M)
        curcell=curcell+1;
    end
    if curcell==M+1
        break;
    end
    poly=[1 t t^2 t^3];
    points(i,:)=poly*curve{1,curcell};
    if i>1
        len=len+norm(points(i,:)-points(i-1,:));
    end
    t=t+dt;
end

end