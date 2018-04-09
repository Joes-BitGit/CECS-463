% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #1 Due: 1/31/2018
clc; clear all; close all; format compact;
fprintf('Lab 1 -- Simple Signals \n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
%Fig 1
N=12; M=4; n=[0:1:2*N-1];
[P,Q]=rat(M/N);
x=sin((2*pi*M*n)/N);
figure(1); hold on; grid on;
axis([0,n(end),-1,1]); xlabel('N'); ylabel('Amp:x(n)')
stem(n,x,'r*');
str=sprintf('Fundamental Period=%d',P); title(str);

%Fig 2
N=12; M=5; n=[0:1:2*N-1];
[P,Q]=rat(M/N);
x=sin((2*pi*M*n)/N);
figure(2); hold on; grid on;
axis([0,n(end),-1,1]); xlabel('N'); ylabel('Amp:x(n)')
stem(n,x,'r*');
str=sprintf('Fundamental Period=%d',P); title(str);

%Fig 3
N=12; M=7; n=[0:1:2*N-1];
[P,Q]=rat(M/N);
x=sin((2*pi*M*n)/N);
figure(3); hold on; grid on;
axis([0,n(end),-1,1]); xlabel('N'); ylabel('Amp:x(n)')
stem(n,x,'r*');
str=sprintf('Fundamental Period=%d',P); title(str);

%Fig 4
N=12; M=10; n=[0:1:2*N-1]; 
[P,Q]=rat(M/N);
x=sin((2*pi*M*n)/N);
figure(4); hold on; grid on;
axis([0,n(end),-1,1]); xlabel('N'); ylabel('Amp:x(n)')
stem(n,x,'r*');
str=sprintf('Fundamental Period=%d',P); title(str);

%Fig 5
N=12; M=15; n=[0:1:2*M-1]; 
[P,Q]=rat(M/N);
x=sin((2*pi*M*n)/N);
figure(5); hold on; grid on;
axis([0,n(end),-1,1]); xlabel('N'); ylabel('Amp:x(n)')
stem(n,x,'r*');
str=sprintf('Fundamental Period=%d',P); title(str);

%% Problem 2
n=0:1:9; k=1;w=2*pi*k/5;
x=sin(w*n);
figure(6); hold on; grid on;
subplot(2,2,1); stem(n,x);
xlabel('n'); ylabel('x[n]');
str=sprintf('x1(n)=sin(%2.2f*n)',w/pi);
title(str);

n=0:1:9; k=2; w=2*pi*k/5;
x=sin(w*n);
figure(6); hold on; grid on;
subplot(2,2,2); stem(n,x);
xlabel('n'); ylabel('x[n]');
str=sprintf('x2(n)=sin(%2.2f*n)',w/pi);
title(str);

n=0:1:9; k=4; w=2*pi*k/5;
x=sin(w*n);
figure(6); hold on; grid on;
subplot(2,2,3); stem(n,x);
xlabel('n'); ylabel('x[n]');
str=sprintf('x3(n)=sin(%2.2f*n)',w/pi);
title(str);

n=0:1:9; k=6; w=2*pi*k/5;
x=sin(w*n);
figure(6); hold on; grid on;
subplot(2,2,4); stem(n,x);
xlabel('n'); ylabel('x[n]');
str=sprintf('x4(n)=sin(%2.2f*n)',w/pi);
title(str);

%% Problem 3
N=6; 
figure(7);
n=linspace(0,71);
T1=N; T2=(2*N/3); L=lcm(T1,T2);
x1=cos(2*pi*n/N)+2*cos(3*pi*n/N);
plot(n,x1);
str=sprintf('x1(n)=cos(2*(pi)*n/%d + 2*cos(3*(pi)*n/%d, Period=%d',T1,T2,L);
title(str);

figure(8);
T3=round(N*pi); T4=round(2*N*pi/3); L=lcm(T3,T4); n=0:12*N-1;
x2=2*cos(2*pi*n/T3)+cos(2*pi*n/T4);
stem(n,x2); xlabel('N'); ylabel('Amplitude');
str=sprintf('x2(n)=cos(2*(pi)*n/%d + 2*cos(2*(pi)*n/%d, Period=%d',T3,T4,L);
title(str);

figure(9);
n=linspace(0,50);
T5=(N); T6=round(4*N/5); L=lcm(T5,T6);
x3=cos(2*pi*n/N)+3*sin((5*pi*n)/2*N);
plot(n,x3);
str=sprintf('x3(n)=cos(2*(pi)*n/%d + 3*sin(5*(pi)*n/%d, Period=%d',T5,T6,L);
title(str);
%% Problem 4
n=0:31; 
x1=sin(pi*n/4).*cos(pi*n/4); figure(10);
stem(n,x1); grid on;
xlabel('n'); ylabel('x1(n)');
str=sprintf('Plot 1'); title(str);

x2=cos((pi*n)/4).^2;
figure(11); hold on; grid on;
stem(n,x2); xlabel('n'); ylabel('x2(n)');
str=sprintf('Plot 2'); title(str);
 
x3=sin((pi*n)/4).*cos((pi*n)/8);
figure(12); hold on; grid on;
stem(n,x3); xlabel('n'); ylabel('x3(n)');
str=sprintf('Plot 3'); title(str);

%% Problem 5
n=0:100; x1=sin(n); x2=sin(n*pi);
figure(13);
subplot(2,1,1); stem(n,x1+x2);
str=sprintf('x(n)=sin(n)+sin(pi*n)'); title(str);

subplot(2,1,2); stem(n,x1.*x2);
str=sprintf('x(n)=sin(n)*sin(pi*n)'); title(str);