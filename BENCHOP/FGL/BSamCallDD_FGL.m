 %Copyright (C) 2015 Magnus Wiktorsson

    %This file is part of BENCHOP.
    %BENCHOP is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %(at your option) any later version.

    %BENCHOP is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.

    %You should have received a copy of the GNU General Public License
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function [U]=BSAmCallDD_FGL(S,K,T,r,sigma,D,alpha)
% function [U]=BSAmCallDD_FGL(S,K,T,r,sigma,D,alpha);
%
% ----------input---------------------------
%
%       S: intial/current stock price
%       K: strike price
%       T: time to maturity
%       r: continously compounded short rate
%   sigma: BM-vol
%       D: relative dividend
%   alpha: divend time as a fraction of time to maturity
%                   i.e. alpha=t_div/T           
%
%   
%  Note that for S,K,T and r that they should either be 1-dim
%  or have matching dimensions
%
%  If t_div>T (i.e. alpha>1) then the price collapses to the standard Black Scholes price
%
% ----------output---------------------------
%
%   U: Black Scholes discrete dividend American call option price
%
% (C) 2014 Magnus Wiktorsson
IS=ones(1,length(S));
if alpha<1
 Sast=binsearch(@(x) x-K-BSeuCallUI_FGL(x*(1-D),K,T*(1-alpha),r,[sigma]),K,2*K,1e-1);
 sast=log(Sast);
 k=log(K);
 Kf=@BSKf;
 load glax50
 [xx1,xx2]=meshgrid(x,x);
 alpha1=6;
 alpha2=6;
 z1=-1i*alpha1+xx1(:);
 z11=-1i*alpha1-xx1(:);
 z2=-1i*alpha2+xx2(:);
 z22=-1i*alpha2-xx2(:);
 ff=@(z1,z2,T,r,s0,sast,sigma,D,alpha)  exp(-1i*z1*k-1i*z2*sast+Kf(-1i+z1+z2,T*alpha,r,s0,sigma)+(1i*z1+1)*log(1-D)+Kf(-1i+z1,T*(1-alpha),r,0,sigma)-log((1i*z2).*(1i*z1+1).*(1i*z1)));
 [w1,w2]=meshgrid(we,we);
 for i=1:length(S)
  s0=log(S(i));
  %f=@(a) (-a(1)*k-a(2)*sast+Kf(-1i-1i*a(1)-1i*a(2),T*alpha,r,s0,sigma)+(a(1)+1)*log(1-D)+Kf(-1i*a(1),T*(1-alpha),r,0*s0,sigma)-log(a(2)*a(1)*(a(1)+1)))*((a(2)>0)&(a(1)>0))+1e34*((a(2)<0)|(a(1)<0));
  %aa(i,:)=fminsearch(f,[1,23],optimset('TolFun',1e-1,'MaxFunEvals',1000));
  %a=[6 6];
  %alpha1=aa(i,1);
  %alpha2=aa(i,2);
  %z1=-1i*alpha1+xx1(:);
  %z11=-1i*alpha1-xx1(:);
  %z2=-1i*alpha2+xx2(:);
  %z22=-1i*alpha2-xx2(:);
  %ff=@(z1,z2,T,r,s0,sast,sigma,D,alpha)  exp(-1i*z1*k-1i*z2*sast+Kf(-1i+z1+z2,T*alpha,r,s0,sigma)+(1i*z1+1)*log(1-D)+Kf(-1i+z1,T*(1-alpha),r,0,sigma)-log((1i*z2).*(1i*z1+1).*(1i*z1)));
  %[w1,w2]=meshgrid(we,we);
  U1(i)=-2*exp(-r*T)*(w1(:).*w2(:))'*real(ff(z1,z2,T,r,s0,sast,sigma,D,alpha)+ff(z1,z22,T,r,s0,sast,sigma,D,alpha))/(2*pi)^2;
 end
 %mean(aa)
 U2=fourieroptspec(Kf,IS,IS,k*IS,sast*IS,log(S),r*IS,T*alpha*IS,sigma,x,we);
 U3=BSeuCallDD_FGL(S,K,T,r,sigma,D,alpha);
 U=U1+U2+U3;
else
  U=BSeuCallDD_FGL(S,K,T,r,sigma,D,alpha); 
end

function [P,ak]=fourieroptspec(Kf,C,p,k,sast,s0,r,Tau,param,x,we)
%function [P,ak]=fourieropt(Kf,C,p,k,s0,r,Tau,param,x,we)
% calculates the optionprice using inverse fouriertransform
% by the Gauss-Laguerre quadrature method, 
% using the optimal straight integration path
% going through the sadlepoint of the integrand
%
%----------- input---------------------------------
%
%         Kf : logchararteristic function for log stockprice
%              use the syntax:
%              [P,ak]=fourieropt(@MertonKf,C,p,k,s0,r,Tau,param,x,we)
%              to use the Merton model and so on
%         C  : row vector containing 1 if call 0 if put option
%         p  : vector containing the p for power options
%              payoff=((S_Tau-K)^+)^p for call i.e. p=1 for an
%              ordinary call option, p=0 for a binary option
%         k  : row vector containg log strike prices
%         s0 : intial/current log stock price 
%         r  : row vector of continously compounded short rate 
%        Tau : row vector containing time to maturity
%       param: model parameters corresponding to the 
%              the chosen logchararteristic function
%              see the help text for the corresponding file
%          x : column vector of evaluation points for the 
%              Gauss-Laguerre quadrature (default using 500 points)
%         we :  column vector of corresponding weights multiplied by exp(x)
%              (default using 500 points)
%          Note that C,p,k,r and Tau should be row vectors with
%          matching dimensions
%
%
%   ----------------output----------------------------
%    P : calculated price
%   ak : vector (matching size(k) ) of offsets used in the 
%         inverse fourier transform
%         kept only for debugging reasons 

% (C) 2006,2014 Magnus Wiktorsson

if nargin==8
 load glax500  
end  
CC=C;
pf=inline('-log(1i*x)','x','C','p');
anintpf=zeros(length(k),2);
if any(C==1)
 anintpf(C==1,:)=ones(sum(C==1),1)*[0,inf];
end
if any(C==0)
 anintpf(C==0,:)=[-inf*ones(sum(C==0),1),-p(C==0)'];
end
[tmp1,anintcf]=Kf(-1i,Tau(1),r(1),s0(1),param);
Tau=Tau(:)';
x=x(:);
we=we(:);
Ik=ones(1,length(k));
Ix=ones(length(x),1);
anint=[max(0.999*max(anintpf(:,1)+sign(anintpf(:,1))*0.01,anintcf(1)-p'),-1e3),min(0.999*min(anintpf(:,2)+sign(anintpf(:,1))*0.01,anintcf(2)-p'),1e3)];
ak=goldsearch_spec(Kf,pf,anint,1e-8,sast,Tau,r(:)',s0(:)',param,C,p);
z=x*Ik-1i*Ix*ak;
P1=exp(-r.*Tau).*(we'*real(exp(Kf(z-1i*Ix*p,Ix*Tau,Ix*r,Ix*s0,param)+pf(z,Ix*C,Ix*p)+(-1i*z.*(Ix*sast)))))/pi;
%ak=goldsearch_spec(Kf,pf,anint,1e-8,sast,Tau,r(:)',s0(:)',param,C,0*p);
%mean(ak)
ak=30.15*Ik;
z=x*Ik-1i*Ix*ak;
P2=exp(-r.*Tau).*(we'*real(exp(Ix*k+Kf(z-1i*Ix*p*0,Ix*Tau,Ix*r,Ix*s0,param)+pf(z,Ix*C,Ix*p)+(-1i*z.*(Ix*sast)))))/pi;
P=(P1-P2);




function [a,d]=goldsearch_spec(Kf,pf,oint,e,k,T,r,s0,par,C,p)
% function [a,d]=goldsearch_spec(Kf,pf,oint,e,k,T,r,s0,par,C,p)
%
% finds the minimum of a convex function by the golden-ratio
% search algorithm  (for internal use in function only)

% (C) 2006 Magnus Wiktorsson
if size(oint,1)==1
 a=oint(1); b=oint(2);
 a=ones(size(k))*a;
 b=ones(size(k))*b;
else
 a=oint(:,1)';
 b=oint(:,2)';
end 
gr=(sqrt(5)-1)/2;
n=ceil(log(e)/log(gr)-log(b(1)-a(1))/log(gr));
m=a+gr*(b-a);
l=a+(1-gr)*(b-a);
fm=(Kf(-1i*(m+p),T,r,s0,par)+pf(-1i*m,C,p)+(-m.*k));
fl=(Kf(-1i*(l+p),T,r,s0,par)+pf(-1i*l,C,p)+(-l.*k));

for i=1:n
 I=(fl>fm); 
 if any(I)
 a(I)=l(I); 
 l(I)=m(I);
 m(I)=a(I)+gr*(b(I)-a(I));
 fl(I)=fm(I);
 fm(I)=(Kf(-1i*(m(I)+p(I)),T(I),r(I),s0(I),par)+pf(-1i*m(I),C(I),p(I))+(-m(I).*k(I)));
 end
 I=~I;
 if any(I)
 b(I)=m(I);
 m(I)=l(I);
 l(I)=a(I)+(1-gr)*(b(I)-a(I));
 fm(I)=fl(I);
 fl(I)=(Kf(-1i*(l(I)+p(I)),T(I),r(I),s0(I),par)+pf(-1i*l(I),C(I),p(I))+(-l(I).*k(I)));
 end
 end
a=(a+b)/2;

function [x]=binsearch(f,a,b,e)
fa=f(a);
fb=f(b);
n=max(ceil((log(b-a)-log(e))/log(2)));
if any((fa.*fb)>0)
   error('Same sign of f at endpoints') % will never happen (I hope)
else   
 for i=1:n
  c=(a+b)/2;
  I=(f(c)<0);
  a(I)=c(I);  
  b(~I)=c(~I);
 end 
end
x=(a+b)/2;


function [K,anint]=BSKf(u,T,r,s0,par)
% log characteristic function for the Black-Scholes model
%
% ------------input-----------------------
%
%     z: evaluation point, Note that we allow z to be complex, if we
%        evaluate at the point -ix (for x real) we evaluate the log
%        moment generating function at the point x. 
%     T: time
%     r: continously compounded short rate 
%    s0: intial/current log stock price 
% param: model parameters
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006 Magnus Wiktorsson

s=par(1); % BM-vol
K=1i*u.*(s0+T.*(r-s.^2/2))-T.*u.^2/2*s.^2;
anint=[-inf,inf];