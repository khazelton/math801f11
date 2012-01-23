function runtorusx(ddd)
su=2*pi/120;
sv=2*pi/90;
u1 = 0:su:(2*pi);
u2 = 0:sv:(2*pi);
for i = 1:10
    ua = 0.015*rand(size(u1));ub=0.015*rand(size(u2));
    u1 = u1 +ua;u1(1)=0;u1(end)=2*pi;
    u2 = u2 +ub;u2(1)=0;u2(end)=2*pi;
    [u,v]=meshgrid(u1,u2);
    uc = 0.01*rand(size(u));
    u = u + uc; v = v + uc;
    r1=2.05; %  radius
    r2=0.36; % Radius
    r3=0.01; % tube radius (going around the torus)
    p1=sqrt(5^2+0.02); % Turns around 1-direction
    p2=sqrt(14^2+0.02); % Twists around 2-direction
    
    x=r1*cos(u*p1) + r2*cos(u*p1).*cos(u*p2) +r3*cos(u*p1).*sin(v);
    y=r1*sin(u*p1) + r2*sin(u*p1).*cos(u*p2) +r3*sin(u*p1).*sin(v);
    z=r2*sin(u*p2)+r3*cos(v);
    %c = zeros(size(z));
    hold on;
    c=z;
    mesh(x,y,z,c,'marker','.','edgecolor','none','facecolor','none','markeredgecolor','b');
end
%colormap colorcube
view(-77,42);
axis equal

figure;
idx2=round(12.2:24:121);
su=2*pi/120;
sv=2*pi/90;
u1 = 0:su:(2*pi);
u2 = 0:sv:(2*pi);
%for j = 1:6
for i = 1:10
    ua = 0.015*rand(size(u1));ub=0.015*rand(size(u2));
    u1 = u1 +ua;u1(1)=0;u1(end)=2*pi;
    u2 = u2 +ub;u2(1)=0;u2(end)=2*pi;
    [u,v]=meshgrid(u1,u2);
    uc = 0.01*rand(size(u));
    u = u + uc; v = v + uc;
    r1=2.05; %  radius
    r2=0.36; % Radius
    r3=0.01; % tube radius (going around the torus)
    p1=sqrt(5^2+0.02); % Turns around 1-direction
    p2=sqrt(14^2+0.02); % Twists around 2-direction
    
    x=r1*cos(u*p1) + r2*cos(u*p1).*cos(u*p2) +r3*cos(u*p1).*sin(v);
    y=r1*sin(u*p1) + r2*sin(u*p1).*cos(u*p2) +r3*sin(u*p1).*sin(v);
    z=r2*sin(u*p2)+r3*cos(v);
    %c = zeros(size(z));
    hold on;
    c=z;
    mesh(x,y,z,c,'marker','.','edgecolor','none','facecolor','none','markeredgecolor','y');
    mesh(x(:,idx2),y(:,idx2),z(:,idx2),c(:,idx2),'marker','p','edgecolor','none','facecolor','none','markeredgecolor','b');
end
mesh(x(:,idx2),y(:,idx2),z(:,idx2),c(:,idx2),'marker','h','edgecolor','none','facecolor','none','markeredgecolor','b');
%end
%colormap colorcube
view(-77,42);
axis equal
torus
