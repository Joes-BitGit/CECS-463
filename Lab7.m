% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #7 Due: 3/5/2018
clc; clear all; close all; format compact;
fprintf('Lab 7 -- The Z-Transform\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
x=[1,0,0,1,-1,3,-2]; nx=[-3:3]; 
y=[1,-2,4,3,-2,1]; ny=[-3:2];
fprintf('Problem 1a: \n');
[x1,n1]=conv_m(x,nx,y,ny)

x2_1=[1,1,1]; n2_1=[0,1,2];
[x2_2,n2_2]=conv_m(x2_1,n2_1,x2_1,n2_1);
fprintf('Problem 1b: \n');
[x2_3,n2_3]=conv_m(x2_2,n2_2,x2_2,n2_2)

[x3_1,n3_1]=sigadd(x,nx,-1*y,ny);
fprintf('Problem 1c: \n');
[x3_2,n3_2]=conv_m(x3_1,n3_1,x3_1,n3_1)

[x4_1,nx4_1]=conv_m(x,nx,x,nx);
[y4_1,ny4_1]=conv_m(y,ny,y,ny);
fprintf('Problem 1d: \n');
[x4_2,n4_2]=conv_m(x4_1,nx4_1,-1*y4_1,ny4_1)

%% Problem 2
clear all;
b=[1,0]; a=[1,-3/5];
[delta,n]=impseq(0,0,2);
x=filter(b,a,delta)
x1=(3/5).^n.*stepseq(0,0,2)
fprintf('Problem 2a Region of Convergence: |z|>(3/5)\n');

b=[1,0]; a=[1,2];
[delta,n]=impseq(0,0,4);
x=filter(b,a,delta)
x2=(-2).^n.*stepseq(-1,-1,3)
fprintf('Problem 2b Region of Convergence: |z|<2\n');

b=[0,0,16]; a=[25,-20,0];
[delta,n]=impseq(0,0,4);
x=filter(b,a,delta)
x3=(0.8).^n.*stepseq(2,0,4)
fprintf('Problem 2c Region of Convergence: |z|>(4/5)\n');

b=[20]; a=[16,-40,25];
[delta,n]=impseq(0,0,4);
x=filter(b,a,delta)
x4=[(n+1).*(5/4).^n].*stepseq(0,0,4)
fprintf('Problem 2d Region of Convergence: |z|>(5/4)\n');

%% Problem 3
a=2/3;
S1=1/(1-(a))
a=-3/5;
S2=1/(1-(a))
a=4/3;
P1=(1-((a)^11))/(1-(a))
a=9/5; n1=5; n2=12;
P2=sums(a,n1,n2)
a1=-5/6; n1=6; n2=9;
a2=4/9; n3=3; n4=8;
P3=sums(a1,n1,n2)+sums(a2,n3,n4)