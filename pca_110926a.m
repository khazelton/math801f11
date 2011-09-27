Mu = [0;0];
Sigma = [1;1];
Npoint = 100;
pdx = ProbDistUnivParam('normal',[Mu(1) Sigma(1)]);
pdy = ProbDistUnivParam('normal',[Mu(2) Sigma(2)]);
xset = pdx.random(1,Npoint);
% xset is a column vector of length 100
% xset values are random samples of a normal distribution with mean 0 and S.D. of 1.
% yset is a column vector in r^100 where y(i) is a random sample from a normal distribution
% with a mean of 0 and S.D. of 1
yset = pdy.random(1,Npoint);
% 'yset is given a level of correlation with xset'
yset = .75*yset + .25*xset
setA = [xset;yset];
save setA setA
A = setA;
figure;scatter(A(1,:),A(2,:),'red')
hold on
axis equal
xlabel('X');ylabel('Y');
grid on
M1 = mean(A,2);
MM = M1*ones(1,size(A,2));
K=A*A'
[V1,D1]=eig(K)
L = V1(:,2)
L=L*ones(1,2)
L(:,2)=-L(:,2)
L=L*6
plot(L(1,:),L(2,:))