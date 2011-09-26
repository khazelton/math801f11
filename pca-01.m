
_________________________
2011-09-26 07:54 GMT-05:00

%-- 2011-09-25 21:47 --%
clc
close all
clear
Mu = [0;0];
Sigma = [1;1];
Npoint = 100;
pdx = ProbDistUnivParam('normal',[Mu(1) Sigma(1)]);
pdy = ProbDistUnivParam('normal',[Mu(2) Sigma(2)]);
xset = pdx.random(1,Npoint);
size(xset)
% xset is a column vector of length 100
% xset values are random samples of a normal distribution with mean 0 and S.D. of 1.
figure;hist(xset,10)
% yset is a column vector in r^100 where y(i) is a random sample from a normal distribution
% with a mean of 0 and S.D. of 1
yset = pdy.random(1,Npoint);
yset = .75*yset + .25*xset
setA = [xset;yset];
save setA setA
size(setA)
figure;scatter(setA(1,:),setA(2,:),'red')
xlabel('X');ylabel('Y');
grid on
hold on
% 'create a working copy of setA called A'
A = setA;
% 'subtract the mean of a row from each element of the row'
M1 = mean(A,2);
MM = M1*ones(1,size(A,2));
size(A)
% 'Let K be the symmetric, positive, semi-definite product of A and A transpose'
K=A*A'
% 'Let V1 be the eigenvectors of K and let D1 be the eigenvalues'
[V1,D1]=eig(K)
% 'so D1(2,2) is the largedt eigenvalue (eigenvalues are on diagonal)
% 'and V1(:,2) is the corresponding eigenvector
L=V1[:,2]
L=V1(:,2)
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
grid on
hold on
figure;scatter(L(1),L(2),'red')
figure;scatter(-L(1),-L(2),'red')
grid on
hold on
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
grid on
hold on
scatter(-L(1),-L(2),'red')
line(-L(1),-L(2))
%-- 2011-09-25 22:18 --%
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
grid on
hold on
clear
Mu = [0;0];
Sigma = [1;1];
Npoint = 100;
pdx = ProbDistUnivParam('normal',[Mu(1) Sigma(1)]);
pdy = ProbDistUnivParam('normal',[Mu(2) Sigma(2)]);
xset = pdx.random(1,Npoint);
size(xset)
% xset is a column vector of length 100
xset(1:20)
% xset values are random samples of a normal distribution with mean 0 and S.D. of 1.
hist(xset,10)
% yset is a column vector in r^100 where y(i) is a random sample from a normal distribution
% with a mean of 0 and S.D. of 1
yset = pdy.random(1,Npoint);
hist(yset,10)
yset = .75*yset + .25*xset
hist(yset,10)
setA = [xset;yset];
save setA setA
size(setA)
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
grid on
hold on
A = setA;
M1 = mean(A,2);
MM = M1*ones(1,size(A,2));
size(A)
K=A*A'
[V1,D1]=eig(K)
L=V1(:,2)
L=-L
line(L(1),L(2))
figure;line(L(1),L(2))
line(L(1),L(2),'Color','r',LineWidth,4)
line(L(1),L(2),'Color','r','LineWidth',4)
figure;line(L(1),L(2),'Color','r','LineWidth',4)
L(1)
figure;plot(L(1),L(2))
figure;plot(L(1),L(2),'-')
figure;plot(0,0,L(1),L(2),'-')
figure;line(0,0,L(1),L(2),'-','LineWidth',4)
figure;line(L(1),L(2),'-','LineWidth',4)
figure;line(L(1),L(2),'LineWidth',4)
L
x=[-.9074,.9074]
y=[-.4202,.4202]
figure;plot(x,y)
xlabel('X');ylabel('Y');
grid on
hold on
figure;scatter(setA(1,:),setA(2,:),'blue')
grid on
hold on
plot(x,y)
x=3*x;y=3*y
x
plot(x,y)
x=2*x
y=2*y
plot(x,y)