% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #8 Due: 3/12/2018
clc; clear all; close all; format compact;
fprintf('Lab 8 -- The Z-Transform Properties\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
b=[0,0,0,3]; a=[1,-1];
[delta,n]=impseq(0,0,3);
x=filter(b,a,delta);
[x1s,ns]=stepseq(3,0,3);

b=[0,0,2,-2]; a=[1,-1];
[delta,n]=impseq(0,0,2);
x=filter(b,a,delta);
[x1i,ni]=impseq(2,0,2);

b=[0,0,2,1]; a=[1,-1];
[delta,n]=impseq(0,0,3);
x=filter(b,a,delta)
x1n=sigadd(2*x1i,ni,3*x1s,ns)
fprintf('Problem 1a Region of Convergence: |z|>1\n');
zplane(-1/2,1)

%% Problem 2
b=[0,0,0.25]; a=[1,-.25,0];
[delta,n]=impseq(0,0,4);
x=filter(b,a,delta)
x=(-0.5).^(n);
h=(-0.5).^(n-2).*stepseq(2,2,6);
x1=sigmult(x,n,h,n)
fprintf('Problem 2a Region of Convergence: |z|>1/4\n');

b=[1,0,0]; a=[1,1,0.25];
[delta,n]=impseq(0,0,4);
x=filter(b,a,delta)
x=(-0.5).^(n+2).*stepseq(-2,-2,2);
h=(-0.5).^(n-2).*stepseq(2,2,6);
x1=conv_m(x,n,h,n)
fprintf('Problem 2b Region of Convergence: |z|>.5\n');
%% Problem 3
b=[1,1,1,1,1,1]; a=[1,2,1];
[delta,n]=impseq(0,0,9);
[d,r]=deconv_m(b,n,a,n)
