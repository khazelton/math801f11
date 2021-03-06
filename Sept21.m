clc
close all
clear

%% sample data points

Mu = [2;3];
Sigma = [2;1];
Npoint = 100;

pdx = ProbDistUnivParam('normal',[Mu(1) Sigma(1)]);
pdy = ProbDistUnivParam('normal',[Mu(2) Sigma(2)]);

xset = pdx.random(1,Npoint);
yset = pdy.random(1,Npoint);
setA = [xset;yset];
save setA setA
figure;scatter(setA(1,:),setA(2,:),'blue')
xlabel('X');ylabel('Y');
grid on
hold on

xset = pdx.random(1,Npoint);
yset = pdy.random(1,Npoint);
setB = [xset;yset];
save setB setB
scatter(setB(1,:),setB(2,:),'red')
xlabel('X');ylabel('Y');
grid on
legend('Set A','Set B');
title('Scatter plot of setA and setB data points');
axis equal

%% SVD computation on setA
A = setA;
M1 = mean(A,2);
MM = M1*ones(1,size(A,2));
Azm = A - MM;
K = Azm*Azm';
[V1,D1] = eig(K);

%% Sampling from setB
Emean = zeros(size(setB,2),1);
Epc = zeros(size(setB,2),1);
Apv = zeros(size(setB,2),1);
for Nsample = 2:size(setB,2)
    SampleB = setB(:,(1:Nsample));
    M2 = mean(SampleB,2);
    MM = M2*ones(1,size(SampleB,2));
    Bzm = SampleB - MM;
    K = Bzm*Bzm';
    [V2,D2] = eig(K);
    Emean(Nsample) = norm(M1 - M2);
    Epc(Nsample) = D1(2,2)/(D1(1,1)+D1(2,2)) - D2(2,2)/(D2(1,1)+D2(2,2));
    Apv(Nsample) = acos(dot(V1(:,2),V2(:,2)))*180/pi;
end
figure;plot(Emean);title('Norm of difference between mean(setA) and mean(sampleB');
xlabel('Number of sample taken from B');
grid on
figure;plot(Epc);title('Difference between first PC of setA and first PC of sampleB');
xlabel('Number of sample taken from B');
grid on
figure;plot(Apv);title('Angle difference (degree) between first PV of setA and first PV of sampleB');
xlabel('Number of sample taken from B');
grid on

%% Randomize Sampling from setB
Nsample = 20;
rng('shuffle') 
[m,n] = size(setB);
p = rand(n,1);
[p1,I] = sort(p);
SampleB = setB(:,I(1:Nsample));
M2 = mean(SampleB,2);
MM = M2*ones(1,size(SampleB,2));
Bzm = SampleB - MM;
K = Bzm*Bzm';
[V2,D2] = eig(K);
norm(M1 - M2)
D1(2,2)/(D1(1,1)+D1(2,2)) - D2(2,2)/(D2(1,1)+D2(2,2))
acos(dot(V1(:,2),V2(:,2)))*180/pi
