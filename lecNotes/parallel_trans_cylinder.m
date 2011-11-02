% parallel transport for a surface that is a cylinder on a curve C
clc;
clear;
\%% defining the curve C, which is a unit cylinder in this case.

\%% defining the symbols.
syms tau;
syms s;
syms t;

\%%% defining an expression that describes s.This equation can be changed so that 
\%%% it meets the desired curve.

\%% I define s1 explictly as it can be used in the future 
\%% computations to show the chain rule of differentiation.
\%% but s1 is nothing but s

s1 = tau+1;

x1 = s;
x2 = sqrt(1 - s^2);
x3 = s+4;
c  = [x1 x2 0];

\%% direction vector

syms c1;
syms c2;
syms c3;
v = [c1 c2 sqrt(1 - c1^2 - c2^2)];
 

\%% parametrizing the curve interms of lambda.

t1 = 2*tau + 3;

gamma =c + t*v;

%diff(gamma,tau);

Ls = diff(c,s)
Lt = v;

N = cross(Ls,Lt)

\%% defining the vector in terms of the basis for TpS 
syms a;
syms b;

eb_symbolic = [a*diff(x1,s)+b*c1 a*diff(x2,s)+b*c2 b*c3]

% Normal
N_alt = [ c3*diff(x2,s)*diff(s1,tau) -c3*diff(x1,s)*diff(s1,tau) (c2*diff(x1,s)-c1*diff(x2,s))*diff(s1,tau) ]

% d(gamma)/dt
V  = [ diff(x1,s)*diff(s1,tau)+c1*diff(t1,tau) diff(x2,s)*diff(s1,tau)+c2*diff(t1,tau) c3*diff(t1,tau)]    

\%% Assuming W as some function of tau. But in actuality W is a parallel vector 
\%% field to the initial curve. Here W is defined explicitly to simplify the
\%% problem.

syms x11;
syms x22;
syms x33;

\%% My assumtion of W, but it could be taken as other expressions also
W1 = x11^2 + 2*x22^2+x33^2;
W2 = 3 * x11^2 + 2*x22^2+x33^2;
W3 = x11^2 + x22^2+x33^2;

\%% Finding the Gradient
grad_W1 = [ diff(W1,x11)*diff(x1,s)*diff(s1,tau) diff(W1,x22)*diff(x2,s)*diff(s1,tau) diff(W1,x33)*diff(x3,s)*diff(s1,tau)];
grad_W2 = [ diff(W2,x11)*diff(x1,s)*diff(s1,tau) diff(W2,x22)*diff(x2,s)*diff(s1,tau) diff(W2,x33)*diff(x3,s)*diff(s1,tau)];
grad_W3 = [ diff(W3,x11)*diff(x1,s)*diff(s1,tau) diff(W3,x22)*diff(x2,s)*diff(s1,tau) diff(W3,x33)*diff(x3,s)*diff(s1,tau)];

\%% Computing the inner product.
inner_W1_V = dot(V,grad_W1);
inner_W2_V = dot(V,grad_W2);
inner_W3_V = dot(V,grad_W3);

\%% Computing Dv.W
DvW = [ inner_W1_V inner_W2_V inner_W3_V]; 

\%% Computing the ODE system . Refer the document for the actual equation.

RHS = dot(DvW,N_alt)*N_alt;

Expr = DvW - RHS

