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

function [U]=BSeuCallGammaII_FGL(S,K,T,r,sigma)
% function [U]=BSeuCallGammaII_FGL(S,K,T,r,sigma);
%
% ----------input---------------------------
%
%       S: intial/current stock price
%       K: strike price
%       T: time to maturity
%       r: continously compounded short rate
%   sigma: BM-vol
%
%   
%  Note that for S,K,T and r that they should either be 1-dim
%  or have matching dimensions
%
% ----------output---------------------------
%
%   U: Black Scholes European call option Gamma
%
% (C) 2014 Magnus Wiktorsson
S=S(:)';
IS=ones(1,length(S));
load glax50
U=fouriergamma(@BSKf,IS,IS,log(K)*IS,log(S),r*IS,T*IS,sigma*IS,x,we);

function [Gamma,ak]=fouriergamma(Kf,C,p,k,s0,r,Tau,param,x,we)
%function [Gamma,ak]=fouriergamma(Kf,C,p,k,s0,r,delta,Tau,param,x,we)
% calculates gamma using inverse 
% fouriertransform by the Gauss-Laguerre quadrature method, 
% using the optimal straight integration path
% going through the sadlepoint of the integrand
%
%----------- input---------------------------------
%
%         Kf : logchararteristic function for log stockprice
%              use the syntax:
%              [P,ak]=fourierdelta(@MertonKf,C,p,k,s0,r,delta,Tau,param,x,we)
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
%    Gamma: calculated gamma of option
%   ak : vector (matching size(k) ) of offsets used in the 
%         inverse fourier transform
%         kept only for debugging reasons 

% (C) 2011 Magnus Wiktorsson

if nargin==8
 load glax500  
end    
pf=inline('betac(p+1,sign(2*C-1).*(i*x)+p.*(C-1),0)+log(1i*x+p)+log(1i*x+p-1)','x','C','p');
if all(p==1)
 pf=inline('zeros(size(x))','x','C','p');
end
if all(p==0)
 pf=inline('-log(sign(2*C-1))','x','C','p');
end
anintpf=zeros(length(k),2);
anintpf(C==1,:)=ones(sum(C==1),1)*[-inf,inf];
anintpf(C==0,:)=[-inf*ones(sum(C==0),1),(-p(C==0).*(1-(p(C==0)==1)))'];
[tmp1,anintcf]=Kf(-1i,Tau(1),r(1),s0(1),param);
   
Tau=Tau(:)';
x=x(:);
we=we(:);
x=x/(param(1)*sqrt(Tau(1)));
we=we/(param(1)*sqrt(Tau(1)));
Ik=ones(1,length(k));
Ix=ones(length(x),1);
anint=[max(0.999*max(anintpf(:,1)+sign(anintpf(:,1))*0.01,anintcf(1)-p'),-1e3),min(0.999*min(anintpf(:,2)+sign(anintpf(:,1))*0.01,anintcf(2)-p'),1e3)];
ak=goldsearch_spec(Kf,pf,anint,5e-1,k,Tau,r(:)',s0(:)',param,C,p); % gamma
%mean(ak)
%ak=-150.9*Ik;
z=x*Ik-1i*Ix*ak;
Gamma = exp(-r.*Tau).*(we'*real(exp(Kf(z-1i*Ix*p,Ix*Tau,Ix*r,Ix*s0,param)-2*Ix*s0+pf(z,Ix*C,Ix*p)+(-1i*z.*(Ix*k)))))/pi;



function a =goldsearch_spec(Kf,pf,oint,e,k,T,r,s0,par,C,p)
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
fm=(Kf(-1i*(m+p),T,r,s0,par)-2*s0+pf(-1i*m,C,p)+(-m.*k));
fl=(Kf(-1i*(l+p),T,r,s0,par)-2*s0+pf(-1i*l,C,p)+(-l.*k));

for i=1:n
 I=(fl>fm); 
 if any(I)
 a(I)=l(I); 
 l(I)=m(I);
 m(I)=a(I)+gr*(b(I)-a(I));
 fl(I)=fm(I);
 fm(I)=(Kf(-1i*(m(I)+p(I)),T(I),r(I),s0(I),par)-2*s0(I)+pf(-1i*m(I),C(I),p(I))+(-m(I).*k(I)));
 end
 I=~I;
 if any(I)
 b(I)=m(I);
 m(I)=l(I);
 l(I)=a(I)+(1-gr)*(b(I)-a(I));
 fm(I)=fl(I);
 fl(I)=(Kf(-1i*(l(I)+p(I)),T(I),r(I),s0(I),par)-2*s0(I)+pf(-1i*l(I),C(I),p(I))+(-l(I).*k(I)));
 end
 end
a=(a+b)/2;

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
