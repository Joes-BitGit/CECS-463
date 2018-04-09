% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #2 Due: 2/5/2018
clc; clear all; close all; format compact;
fprintf('Lab 2 -- Discrete Time Signals \n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
figure(1);
n=-3:7;
x=2*impseq(0,-3,7)+impseq(2,-3,7)+-1*impseq(3,-3,7)+3*impseq(4,-3,7);
stem(n,x); grid on;
xlabel('n'); ylabel('d(n)')
str=sprintf('Discrete Time Signal d[n]'); title(str);
    
%% Problem 2
figure(2);
n=-10:10;
x=n;
[y1,m1]=sigshift(x,n,2);
subplot(2,2,1);stem(m1,y1);
title('y1[n]=x[n-2]');
xlabel('Index'); ylabel('y1[n]');

[y2,m2]=sigshift(x,n,-1);
subplot(2,2,2);stem(m2,y2);
title('y2[n]=x[n+1]');
xlabel('Index'); ylabel('y2[n]');

[y3,m3]=sigfold(x,n);
subplot(2,2,3);stem(m3,y3);
title('y3[n]=x[-n]');
xlabel('Index'); ylabel('y3[n]');

[y4,m4]=sigfold(x,n);
[y4,m4]=sigshift(y4,m4,-1);
subplot(2,2,4);stem(m4,y4);
title('y4[n]=x[-n+1]');
xlabel('Index'); ylabel('y4[n]');
