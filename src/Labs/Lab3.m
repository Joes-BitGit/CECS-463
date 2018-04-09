% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #3 Due: 2/7/2018
clc; clear all; close all; format compact;
fprintf('Lab 3 -- Linear Shift Invariant (LSI) Systems\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1a
% y(n)=T[x(n)]=sin((pi/2)*x(n))
A=2; B=-1; %Arbitrary Constants
x1=impseq(0,0,10); x2=2*impseq(0,0,10); %Arbitrary Sequences
n=0:10;
figure(1); hold on; grid on;
subplot(2,1,1); stem(n,x1,'b');
xlabel('n'); ylabel('x1(n)');
subplot(2,1,2); stem(n,x2,'r');
xlabel('n'); ylabel('x2(n)');
%Linearity: y[Ax1[n]+Bx2[n]]=Ay[x1[n]]+By[x2[n]]
[ylhs,nlhs]=sigadd(A*(pi/2)*x1,n,B*(pi/2)*x2,n);
ylhs=sin(ylhs);
x1f=sin((pi/2)*x1); nx1=n; x2f=sin((pi/2)*x2); nx2=n;
[yrhs,nrhs]=sigadd(A*x1f,nx1,B*x2f,nx2);
figure(2);
subplot(2,1,1); grid on; stem(nlhs,ylhs,'b');
title('LHS = T[Ax1[n]+Bx2[n]]'); xlabel('n'); ylabel('LHS');
subplot(2,1,2); grid on; stem(nrhs,yrhs,'r');
title('RHS = AT[x1[n]]+BT[x2[n]]'); xlabel('n'); ylabel('RHS');
if(sum(abs(yrhs-ylhs))>10*eps)
    fprintf('y(n)=sin((pi/2)*x(n)) is NOT LINEAR\n');
else
    fprintf('y(n)=sin((pi/2)*x(n)) is LINEAR\n');
end

%% Problem 1b
% y(n)=T[x(n)]=n*x(n)
n=0:10; x=impseq(0,0,10); k=-1;
ylhs=n.*x; [ylhs,nlhs]=sigshift(ylhs,n,k);
[yrhs,nrhs]=sigshift(x,n,k); yrhs=nrhs.*yrhs;
figure(3); grid on;
subplot(2,1,1); stem(nlhs,ylhs,'b');
title('LHS = S[T[x(n)]]'); xlabel('n'); ylabel('LHS');
subplot(2,1,2); stem(nrhs,yrhs,'r');
title('RHS = T[S[x(n)]]'); xlabel('n'); ylabel('RHS');
[diff,ndiff]=sigadd(ylhs,nlhs,-yrhs,nrhs);
if(sum(abs(diff))>10*eps)
    fprintf('y(n)=T[x(n)]=n*x(n) is NOT SHIFT INVARIANT\n');
else
    fprintf('y(n)=T[x(n)]=n*x(n) is SHIFT INVARIANT\n');
end

%% Problem 1c
% y(n)=x(n)+x(n+1)
n=-10:10; xn=-5:9; yn=-6:9;
x=stepseq(0,-5,9);
y=stepseq(0,-6,9)+stepseq(-1,-6,9);
figure(4);
axis([-10,10, min(y)-1,max(y)+1]);
subplot(2,1,1); stem(x);
title('x(n)');
subplot(2,1,2); stem(y);
title('y(n)');
