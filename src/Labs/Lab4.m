% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #5 Due: 2/14/2018
clc; clear all; close all; format compact;
fprintf('Lab 5 -- The Fourier Series\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
T=1; %Fundamental Period
F=1/T;
dT=0.01; %step size
t = -1.1:0.01:1.1; %declare the range of t over the fundamental period 0<t<1
x=.5*sign(sin(2*pi*F*t)); %square wave
%sq=.5*(square(2*pi*F*t,50)); %square wave
figure();
plot(t,x,'b');
hold on; grid on;
axis([-1.1, 1.1, -1.1, 1.1]); title('Periodic Square Wave'); 
xlabel('time t');ylabel('f(t)');
c0=1/T*sum(x*dT);
N=1;
for n=1:N
    cn(n)=dT/T*sum(x.*exp(-1j*2*pi*n*t/T));
end
c_n=conj(cn);

Wn=exp(1j*2*pi/T * t'*n); % T is period, t is time vector, n is coefficient index vector. 
W_n=conj(Wn); % Simply the conjugate of Wn, & c0 is average value of function 
gN=(c0+Wn*cn).'+(W_n*c_n).'; % Careful not to conjugate transpose the row vector cn 

plot(t,gN,'r');