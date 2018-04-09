% MATLAB HEADER
% Joseph Almeida
% CECS 463 SOC II Sp18
% Assignment #10 Due: 3/21/2018
clc; clear all; close all; format compact;
fprintf('Lab 10 -- The Z-transform\n');
str=datestr(now); fprintf('Matlab Time Stamp: %s\n',str);

%% Problem 1
fprintf('Problem 1\n');
n=0:20;
x1=(2/3).^n.*(n>=0);
x2=-(3/4).^n.*(n>=4);
x3=-(1/2).^n.*(n==2);
x=x1+x2+x3;
%outVector(x,n,'x');
figure(); subplot(2,1,1); stem(n,x); grid on;
title('Problem 1: x(n)=[ (2/3)^{n}-(3/4)^{n}-(1/2)^{n} ] u(n)');
xlabel('n'); ylabel('x(n)');
%n=0:20; b=[1]; a=[1,-2/3]; y1=filter(b,a,n==0);
%n=0:20; b=[0,0,0,0,-(3/4)^4]; a=[1,-3/4]; y2=filter(b,a,n==0);
%n=0:20; b=[0,0,-(1/2)^2]; a=1; y3=filter(b,a,n==0);
%outVector(y1+y2+y3,n,'y');
%figure(); stem(n,y1+y2+y3);
a=conv([1,-2/3],[1,-3/4]);
b1=conv(a,[0,0,(1/2)^2]);
b2=conv([1,-2/3],[0,0,0,0,(3/4)^4]);
b3=[1,-3/4];
[b,nb]=sigadd(b1,[0:length(b1)-1],b2,[0:length(b2)-1]);
[b,nb]=sigadd(b3,[0:length(b3)-1],-b,nb);
y=filter(b,a,n==0);
subplot(2,1,2); stem(n,y); grid on;
title('Problem 1: x(n)=Z^{-1}[X(z)]');
xlabel('n'); ylabel('x(n)');
outVector(y,n,'y');

%% Problem 2
n=0:20; w0=pi/3; s=sin(w0); c=cos(w0); p=(6/5)^-1;
x=n.*sin(w0*n)+p.^n.*(n>=2);
%outVector(x,n,'x');
fprintf('  X1(z) = -z(d/dz)[sin(w0)z^-1/(1-2cos(w0)z^-1+z^-2)] = sin(w0)(z^-1+z^-3)/(1-2cos(w0)z^-1+z^-2))^2\n');
fprintf('  X2(z)= (5/6)^2 z^-2/(1-(5/6)z^-1)\n');
bx1=[0,s,0,-s]; nbx1=0:3;
ax1=conv([1,-2*c,1],[1,-2*c,1]); nax1=0:4;
bx2=[0,0,p^2]; nbx2=0:2;
ax2=[1,-p]; nax2=0:1;
[n1,nn1]=conv_m(bx1,nbx1,ax2,nax2); [n2,nn2]=conv_m(bx2,nbx2,ax1,nax1);
[b,m]=sigadd(n1,nn1,n2,nn2);
[a,n]=conv_m(ax1,nax1,ax2,nax2);
n=0:20; y=filter(b,a,n==0);
fprintf('  X(z) = X1(z)+X2(z) = b(z)/a(z)\n');
fprintf('  [');
for k=1:length(b) fprintf('(%6.4f)z^(-%d)+',b(k),k-1); end;
fprintf('\b]\n  =  -------------------------------------------------------------------------------------------------------------\n');
fprintf('  [');
for k=1:length(a) fprintf('(%6.4f)z^(-%d)+',a(k),k-1); end;
fprintf('\b]\n\n');
fprintf('  Error: x1(n)-Z^{-1}[X1(z)] = %6.4f\n',sum(abs(x-y)));
r=roots(a), m=1:length(r); k=min(m(abs(r)==max(abs(r))));
fprintf('ROC: |z|>%4.2f\n',abs(r(k)));
%outVector(b,m,'b',6);
%outVector(a,n,'a',6);
%outVector(y,n,'y');
figure(); stem(n,y);

%% Problem 3
b1=[0,0,0,1]; m1=[0:3]; b2=[0,0,0,0,1]; m2=[0:4]; f=0.75
a1=[1,-f]; n1=0:1; a2=[-f, 1]; n2=0:1;
[a,n]=conv_m(a1,n1,a2,n2);
[b11,nb1]=conv_m(a1,n1,b2,m2);
[b12,nb2]=conv_m(b1,m1,a2,n2);
[b,m]=sigadd(b11,nb1,b12,nb2);
%outVector(b,m,'b');
%outVector(a,n,'a');
fprintf('  Y(z)=Z[x(n-3)+x(3-n)]=Z[x(n-3)+x(-(n-3))] = z^-3 X(z)+ z^-3 X(1/z)\n');
fprintf('  = z^-3[ 1/(1-0.75z^-1) + 1/(1-0.75z)] with ROC 3/4<|z|< 4/3\n');
fprintf('  ');
for k=1:length(b)
if(b(k)~=0) fprintf('(%6.4f)z^(-%d)+',-b(k)/f,k-1); end;
end;
fprintf('\b \n');
fprintf('  =  -----------------------------------------------\n')
fprintf('  ');
for k=1:length(a)
fprintf('(%6.4f)z^(-%d)+',-a(k)/f,k-1);
end;
fprintf('\b \n');
%fprintf('  Error=%6.4f\n',sum(abs(x-y)));
n=-20:20;
b=[1]; a=[1,-0.75]; x3=filter(b,a,n==3); %Shift by 3->impulse at n=3
b=[0,1]; a=[-0.75,1]; x_3=antifilter(b,a,n==3);  %Shift by 3->impulse at n=3
[y,ny]=sigadd(x3,n,x_3,n);
figure(); subplot(2,1,1); stem(ny,y); axis([-20,20,0,3]); grid on;
title('Problem 3: y(n) = Z^{-1}[ Z[x(n-3)+x(3-n)] ] where x(n)=(3/4)^{n}u(n)');
xlabel('n');ylabel('y(n)');
N=20; n=-N:N; f=0.75;
[u3,n3]=sigshift(n>=0,n,3);
[uf,nf]=sigfold(n>=0,n); [u_3,n_3]=sigshift(uf,nf,3);
%figure(); stem(n3,u3)
%figure(); stem(n_3,u_3)
[x1,nx1]=sigmult(f.^(n-3),n,u3,n3);  %outVector(x1,nx1,'x1');
[x2,nx2]=sigmult(f.^(3-n),n,u_3,n_3);  %outVector(x2,nx2,'x2');
[x,nx]=sigadd(x1,nx1,x2,nx2);  %outVector(x,nx,'x');
subplot(2,1,2); stem(nx,x); axis([-N,N,0,3]); grid on;
title('Problem 3: y(n) = x(n-3) + x(3-n) where x(n)=(3/4)^{n}u(n)');
xlabel('n');ylabel('y(n)');
%subplot(3,1,1); stem(nx1,x1); axis([-N,N,0,4]);
%subplot(3,1,2);stem(nx2,x2); axis([-N,N,0,4]);
%subplot(3,1,3); stem(nx,x); axis([-N,N,0,4]);

%% Problem 4
fprintf('x(n)=(1/2)^n u(n)\n');
fprintf('X1(z) = zX(1/z) so x1(n)= Z^-1[X(1/z)] for n->n+1\n');
fprintf('x1(n)= (1/2)^-n u(-n) for n->n+1 so x1(n)= 2^(n+1) u(-n-1)\n');

%% Problem 5
clear all;
fprintf('Problem 5\n');
fprintf(' x(z)/y(z)=p(z)+r(z)/y(z)\n');
n=-2:3; x=[1,0,1/2,-2,0,-4/9]; [x1,nx1]=conv_m(x,n,x,n);
n=-1:2; x=[1,4/3,3,-8/7]; [x2,nx2]=conv_m(x,n,x,n); [x2,nx2]=conv_m(x2,nx2,x,n);
[b,n]=sigadd(x1,nx1,x2,nx2);
m=-2:4; a=[1,0,1/6,3,3/4,0,3/8];
[p,np,r,nr]=deconv_m(b,n,a,m);
outVector(b,-n,'b'); outVector(a,-m,'a');
outVector(p,-np,'p'); %outVector(np,np,'np');
outVector(r,-nr,'r'); %outVector(nr,nr,'nr');

%% Problem 6
n=0:10;
b=[1,-1,-4,4]; a=[1,-11/4,13/8,-1/4];
[r,p,C]=outRPC(b,a); y=outRPCN(r,p,C,n);
outVector(y,n,'y');
figure(); subplot(2,1,1); stem(n,y); grid on;
title('Problem 6(a): Output x(n)  ROC: |z|>2');xlabel('n');ylabel('x(n)');
fprintf('Problem 6(b)\n');
n=0:10;
b=[1,1,-4,4]; a=[1,-11/4,13/8,-1/4];
[r,p,C]=outputRPC(b,a);
x1=antifilter(-r(1),[1,-p(1)],n); outVector(x1,n,'x1(ls)');
x2=yRPCN(r(2:3),p(2:3),C,n); outVector(x2,n,'x2(rs)');
x=x1+x2;
subplot(2,1,2); stem(n,x); grid on;
title('Problem 6(b): Output x(n) ROC: 0.5<|z|<2');xlabel('n');ylabel('x(n)');

%% Problem 7
n=0:20; b=[1,1]; a=[1,-0.5];
[r,p,C]=outRPC(b,a); h=outRPCN(r,p,C,n);
figure(); stem(n,h); grid on; axis([0,20,0,2]);
title('Problem 7: Impulse Response h(n)'); xlabel('n'); ylabel('h(n)');
fprintf('h(n)=Z^-1[1/(1-0.5z^-1)]+Z^-1[z^-1/(1-0.5z^-1)]\n  =(0.5)^n u(n)+(0.5)^(n-1) u(n-1)\n');
fprintf('y(n)=x(n)+x(n-1)+0.5y(n-1)\n');

%% Problem 8
n=0:25; b=[1,0.5]; a=[1,0.5,-0.25];
[r,p,C]=residuez(b,a);
fprintf('h(n)= ');
for k=1:length(r)
fprintf('(%6.4f)(%6.4f)^n+',r(k),p(k));
end
fprintf('\b\n');
fprintf('H(z)= (1+0.5z^-1)/(1+0.5z^-1-0.25z^-2)\n');
b=2*b; a=conv(a,[1,-0.9]);
[r,p,C]=outputRPC(b,a);
fprintf('y(n)= ');
for k=1:length(r)
fprintf('(%6.4f)(%6.4f)^n+',r(k),p(k));
end
fprintf('\b\n');
y=yRPCN(r,p,C,n); figure(); stem(n,y); axis([0,20,0,2.5]); grid on;
title('Problem 8: Output Response y(n) to input x(n)=2(0.9)^{n} u(n)'); xlabel('n'); 
ylabel('y(n)');

%% Problem 9
n=0:50; f=1/2; x=(f).^n.*(n>=0);
a=[1,-0.4,-0.45]; b=[0.45,0.4,-1]; bi=1; ai=[1,-f];
Y=[0,3]; X=[2,2]; xic=filtic(b,a,Y,X)
bf=b+conv(xic,ai)
af=conv(ai,a)
[r,p,C]=outputRPC(bf,af);
fprintf('h(n)= ');
for k=1:length(r)
fprintf('(%6.4f)(%6.4f)^n+',r(k),p(k));
end
fprintf('\b\n');
ytr=outRPCN(r,p,C,n); yss=0*n; ytot=yss+ytr;
figure(); hold on;
subplot(2,1,1); hold on; stem(n,yss,'b'); stem(n,ytr,'r'); grid on; legend('yss','ytr');
title('Problem 9: y_{ss} and y_{tr} solutions'); xlabel('n'); ylabel('Amplitude'); 
subplot(2,1,2); stem(n,ytot,'k'); grid on;
title('Problem 9: y_{total} =y_{ss}+y_{tr} solutions'); xlabel('n'); ylabel('ytot(n)');
        

