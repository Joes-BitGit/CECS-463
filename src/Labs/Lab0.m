% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #0 Due: 1/31/2018
clc; clear all; close all; format compact;
fprintf('Lab 0 -- The Complex Plane \n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);
%%Problem 1
%Complex Operations
%a
x=1-2j; y=-1+1j;
a=(x+y); b=(x-y); c=(x*y); d=(x/y); p=x^y;

%b
figure(1); hold on; axis([-6,6,-6,6]); grid on;
ar=real(a); ai=imag(a);
br=real(b); bi=imag(b);
cr=real(c); ci=imag(c);
dr=real(d); di=imag(d);
pr=real(p); pi=imag(p);
plot(ar,ai,'g*');
plot(br,bi,'k*');
plot(cr,ci,'y*');
plot(dr,di,'c*');
plot(pr,pi,'m*');
title('Results of Complex Operations');
xlabel('REAL'); ylabel('IMAGINARY');
legend('x+y','x-y','x*y','x/y','x^+y');
