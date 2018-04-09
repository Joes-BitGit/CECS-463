% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #6 Due: 2/26/2018
clc; clear all; close all; format compact;
fprintf('Lab 6 -- The Discrete Time Fourier Transform(DTFT)\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
N=1000; n=0:100; h=(0.85).^n; k= 0:N; w=2*pi*k/N;
H=h*exp(-1j*(n'*w)); Hm=abs(H); Ha=angle(H);
figure(1);
subplot(2,1,1); plot(w/pi,Hm); grid on;
title('The DTFT of h(n)=0.85^n'); xlabel('\omega/\pi'); ylabel('|H(\omega)|');

subplot(2,1,2); plot(w/pi,Ha*180/pi); grid on;
title('The DTFT of h(n)=0.85^n'); xlabel('\omega/\pi'); ylabel('\angleH(\omega)');
fprintf('Problem 1: See Figure 1\n');

%% Problem 2
w0=0.2*pi; x=cos(w0*n); 
figure(2);

subplot(2,1,1); stem(n,x); grid on;
title('Input x(n)'); xlabel('n'); ylabel('x(n)');

y=conv(x,h); y=y((1:length(x))); %Truncates the output
subplot(2,1,2); stem(n,y); grid on;
title('Output y(n)'); xlabel('n');ylabel('x(n)');
fprintf('Problem 2: y(n)= |H(w==w0)| cos( w0 (n+ Ha(w==w0)/w0 )\n');

gain=Hm(w==w0);
fprintf(' Gain = %4.2f \n',gain);

phi=Ha(w==w0)/w0;
fprintf(' Phase shift = %4.2f samples\n',phi);

%% Problem 3
n=0:20; x=sin(0.25*pi*n); k=3; [xs,ns]=sigshift(x,n,k);
N=1000; w=2*pi*[0:N]/N;

Xs=xs*exp(-1j*ns'*w); Xsm=abs(Xs); Xsa= angle(Xs);
X=(x*exp(-1j*n'*w)).*exp(-1j*w*k); Xm=abs(X); Xa=angle(X);

MaxError=max(abs(X-Xs));
fprintf('Problem 3: Max Error(LHS-RHS = %6.4f) (see Figure 3)\n',MaxError);
figure(3);

subplot(2,2,1); plot(w/pi,Xsm); grid on; 
title('LHS: DTFT[x(n-k)]'); xlabel('\omega/\pi'); ylabel('|X(\omega)|');

subplot(2,2,3); plot(w/pi,Xsa*180/pi); grid on;
title('The DTFT[x(n-k)]'); xlabel('\omega/\pi'); ylabel('\angleX(\omega)');

subplot(2,2,2); plot(w/pi,Xm); grid on; 
title('RHS: DTFT[x(n)] e^{-jwk}'); xlabel('\omega/\pi'); ylabel('|X(\omega)|');

subplot(2,2,4); plot(w/pi,Xa*180/pi); grid on;
title('The DTFT[x(n)]'); xlabel('\omega/\pi'); ylabel('\angleX(\omega)');
%% Problem 4
clear all; 
N=1000; w=[-N:N]*pi/N; a=[1,0.8,-0.95]; na=0:2; b=[0.5]; nb=0:0;
H=(a*exp(-1j*na'*w))./(b*exp(-1j*nb'*w)); Hm=abs(H); Ha=angle(H);
figure(4);

subplot(2,1,1); plot(w/pi,Hm); grid on;
title('Problem 4: The DTFT H(\omega)'); xlabel('\omega/\pi'); ylabel('|H(\omega)|');

subplot(2,1,2); plot(w/pi,Ha*180/pi); grid on;
title('Problem 4: The DTFT H(\omega)'); xlabel('\omega/\pi'); ylabel('\angleH(\omega)');

fprintf('Problem 4: See Figure 4\n');
%% Problem 5
w0=0.05*pi;
fprintf('Problem 5: y(n)=|H(w0)|cos(w0(n+Ha(w0)/w0 for w0=%4.2f\x03c0 rad/sample\n', w0/pi);

gain=Hm(w==w0);
fprintf(' Gain = %4.2f \n',gain);

phi=Ha(w==w0)/w0;
fprintf(' Phase shift = %4.2f samples\n',phi);