%LECTURE #1 Program
%Example of Rotating Phasor in Complex Plane

clear all; close all; format compact; clc;
str = datestr(now); fprintf('Matlab Time Stamp: %s\n', str);

%%Part 1
%Generate vector to represent phasor in complex plane

figure();   %creates a new figure window, and returns its handle
hold on;    %holds the current plot and all axis properties, including 
            %the current color and linestyle, so that subsequent graphing commands
            %add to the existing graph without resetting the color and linestyle.
x=2+2j; xr=real(x); xi=imag(x); xm=round(abs(x))+1; realx(1)=xr;
axis([-xm,xm,-xm,xm]);%([XMIN XMAX YMIN YMAX]) sets scaling for the x- and y-axes
                      % on the current plot
grid on; %adds major grid lines to the current axes.
arrow([0,0],[xr,xi],16,'BASEANGLE',60,'Color', 'k'); %arrow([1 2 3],[0 0 0],36,'BaseAngle',60) 
                                                     %creates an arrow from (1,2,3) to
                                                     %the origin, with an arrowhead of length 36 pixels 
                                                     %and 60-degree base angle.
arrow([0,0],[xr,0],16,'BASEANGLE',60,'Color', 'r');
title('Phasor x(t)=cos(2/PIft+45^o)'); xlabel('REAL');ylabel('IMAG');
T=0.1; f=1/T; dT=T/10; fs=1/dT; colors=['k','b','r','g','y','m','c'];
pause(0.2);
for n=1:50
    arrow([0,0],[xr,xi],16,'BASEANGLE',60,'Color', 'w');%Erase old vector
    arrow([0,0],[xr,0],16,'BASEANGLE',60,'Color', 'w');%Erase old vector
    x=x*exp(1j*2*pi*f*dT); xr=real(x); xi=imag(x); realx(n+1)=xr;
    arrow([0,0],[xr,xi],16,'BASEANGLE',60,'Color', 'k');%Plot new vector
    arrow([0,0],[xr,0],16,'BASEANGLE',60,'Color', 'r');%Plot new vector

    pause(0.25);
end

%%Part 2
%Plot the samples x(nT) and the continous waveform x(t)
figure(); hold on;
plot([0:n]*dT,realx,'ro-'); grid on;
title('Plot of real componenet of phasor');xlabel('time(sec)');ylabel('Amplitude');
%t=0:0.0001:n*dT; xa=abs(x)*cos(2*pi*f*t+angle(x)); plot(t,xa,'b');


