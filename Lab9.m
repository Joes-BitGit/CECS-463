% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #9 Due: 3/14/2018
clc; clear all; close all; format compact;
fprintf('Lab 9 -- The System Function H(z)\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
b=[1,-1,-4,4];
a=[1,-2.75,1.625,-.25];
[r,p,C]=residuez(b,a)

b=[0,0,1]; 
a=[1,2,1.25,.25];
[r,p,C]=residuez(b,a)

%% Problem 2
fprintf('Problem 2a(i) System Representation: 5/(1-(1/4)z^-1)\n');
fprintf('Problem 2a(ii) Difference Equation y(n): y(n)=1/4y(n-1)+5x(n)\n');
b=[5]; a=[1,-1/4];
figure();
zplane(b,a)

fprintf('Problem 2b(i) System Representation H(z): (1-z^-1+2z^-2)/(1-z^-1+z^-2-z^-3)\n');
fprintf('Problem 2b(ii) Difference Equation y(n): y(n)=y(n-1)-y(n-1)+y(n-3)+x(n)-x(n-1)+2x(n-2)\n');
b=[1,-1,2]; 
a=[1,-1,1,-1];
figure();
zplane(b,a)

%% Problem 3
b=[1,-1]; a=[1,0,-0.81];
[r,p,C]=residuez(b,a)
%y=r'*p^n

        

