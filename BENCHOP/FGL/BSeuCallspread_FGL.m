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

function [U]=BSeuCallspread_FGL(S,T,r,sig1,sig2,rho)
% function [U]=BSeuCallspread_FGL(S,K,T,r,sig1,sig2,rho)
%
% ----------input---------------------------
%
%       S: [intial/current stock price 1 ; intial/current stock price 2]
%          (2xN -matrix)
%       K: Strike price
%       T: time to maturity
%       r: continously compounded short rate 
%     sig1: vol stock 1
%   sigma2: vol stock 2
%      rho: correlation between BM:s
%
%  Note that for S1,S2 and T that they should either be 1-dim
%  or have matching dimensions
%
% ----------output---------------------------
%
%   U: Black Scholes European spread option price
%
% (C) 2014 Magnus Wiktorsson
S=S';
S1=S(1,:);
S2=S(2,:);
IS=ones(1,length(S1));
load glax50
U=fourieropt2d(@BS2dKf,IS,IS,log(S1),log(S2),r*IS,T*IS,[sig1,sig2,rho],x,we);

function [P,ak]=fourieropt2d(Kf,C,p,s01,s02,r,Tau,param,x,we)
%function [P,ak]=fourieropt2d(Kf,C,p,s01,s02,r,Tau,param,x,we)
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
%         s01 : intial/current log stock1 price 
%         s02 : intial/current log stock1 price 
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
%C=ones(size(C));
pf=inline('betac(p+1,sign(2*C-1).*(i*x)+p.*(C-1),0)','x','C','p');
if all(p==1)
 pf=inline('-log(1i*x.*(1i*x+1))','x','C','p');
end
if all(p==0)
 pf=inline('-log(1i*x.*sign(2*C-1))','x','C','p');
end
anintpf=zeros(length(s01),2);
if any(C==1)
 anintpf(C==1,:)=ones(sum(C==1),1)*[0,inf];
end
if any(C==0)
 anintpf(C==0,:)=[-inf*ones(sum(C==0),1),-p(C==0)'];
end
[tmp1,anintcf]=Kf(-1i,-1i,Tau(1),r(1),s01(1),s02(1),param);
Tau=Tau(:)';
x=x(:);
we=we(:);
Ik=ones(1,length(s01));
Ix=ones(length(x),1);
anint=[max(0.999*max(anintpf(:,1)+sign(anintpf(:,1))*0.01,anintcf(1)-p'),-1e3),min(0.999*min(anintpf(:,2)+sign(anintpf(:,1))*0.01,anintcf(2)-p'),1e3)];
%ak=goldsearch_spec2d(Kf,pf,anint,1e-8,Tau,r(:)',s01(:)',s02(:)',param,C,p);
%mean(ak)
ak=12*Ik;
z=x*Ik-1i*Ix*ak;
P=exp(-r.*Tau).*(we'*real(exp(Kf(z-1i*Ix*p,-z,Ix*Tau,Ix*r,Ix*s01,Ix*s02,param)+pf(z,Ix*C,Ix*p))))/pi;
%P(~CC)=P(~CC)-(exp(s0(~CC))-exp(k(~CC)-r(~CC).*Tau(~CC)));

function [a,d]=goldsearch_spec2d(Kf,pf,oint,e,T,r,s01,s02,par,C,p)
% function [a,d]=goldsearch_spec2d(Kf,pf,oint,e,T,r,s01,s02,par,C,p)
%
% finds the minimum of a convex function by the golden-ratio
% search algorithm  (for internal use in function only)

% (C) 2006, 2014 Magnus Wiktorsson
if size(oint,1)==1
 a=oint(1); b=oint(2);
 a=ones(size(s01))*a;
 b=ones(size(s01))*b;
else
 a=oint(:,1)';
 b=oint(:,2)';
end 
gr=(sqrt(5)-1)/2;
n=ceil(log(e)/log(gr)-log(b(1)-a(1))/log(gr));
m=a+gr*(b-a);
l=a+(1-gr)*(b-a);
fm=(Kf(-1i*(m+p),1i*m,T,r,s01,s02,par)+pf(-1i*m,C,p));
fl=(Kf(-1i*(l+p),1i*l,T,r,s01,s02,par)+pf(-1i*l,C,p));

for i=1:n
 I=(fl>fm); 
 if any(I)
 a(I)=l(I); 
 l(I)=m(I);
 m(I)=a(I)+gr*(b(I)-a(I));
 fl(I)=fm(I);
 fm(I)=(Kf(-1i*(m(I)+p(I)),1i*m(I),T(I),r(I),s01(I),s02(I),par)+pf(-1i*m(I),C(I),p(I)));
 end
 I=~I;
 if any(I)
 b(I)=m(I);
 m(I)=l(I);
 l(I)=a(I)+(1-gr)*(b(I)-a(I));
 fm(I)=fl(I);
 fl(I)=(Kf(-1i*(l(I)+p(I)),1i*l(I),T(I),r(I),s01(I),s02(I),par)+pf(-1i*l(I),C(I),p(I)));
 end
 end
a=(a+b)/2;

function [K,anint]=BS2dKf(u1,u2,T,r,s01,s02,par)
% log characteristic function for the 2d-Black-Scholes model
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
% (C) 2006, 2014 Magnus Wiktorsson

sig1=par(1); % S1-vol
sig2=par(2); % S2-vol
rho=par(3); % Correlation
K=1i*u1.*(s01+T.*(r-sig1.^2/2))-T.*u1.^2/2*sig1.^2+1i*u2.*(s02+T.*(r-sig1.^2/2))-T.*u2.^2/2*sig2.^2-T.*u1.*u2*sig1*sig2*rho;
anint=[-inf,inf];
