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

function [U]=HSTeuCall_FGL(S,K,T,r,V,kap,th,sig,rho)
%function [U]=HSTeuCall_FGL(S,K,T,r,V,kap,th,sig,rho)
% Calculates the Heston European call optionprice using inverse fouriertransform
% by the Gauss-Laguerre quadrature method, 
% using the optimal straight integration path
% going through the saddlepoint of the integrand
%
%----------- input---------------------------------
%           S: intial/current stock price 
%           K: row vector containg strike prices
%           T: row vector containing time to maturity
%           r: row vector of continously compounded short rate 
%           V: Initial vol.
%         kap: Mean reversion
%          th: Mean level
%         sig: Vol diffusion.
%         rho: Correlation
%
%          Note that K,S, r and T should be row vectors with
%          matching dimensions
%
%
%   ----------------output----------------------------
%    U : calculated call price
% (C) 2014 Magnus Wiktorsson
S=S(:)';
IS=ones(1,length(S));
load glax50
U=fourieropt(@HestonKf,IS,IS,log(K)*IS,log(S),r*IS,T*IS,[V,kap,th,sig,rho],x,we);

function [P,ak]=fourieropt(Kf,C,p,k,s0,r,Tau,param,x,we)
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
%C=ones(size(C));
pf=inline('betac(p+1,sign(2*C-1).*(i*x)+p.*(C-1),0)','x','C','p');
if all(p==1)
 pf=inline('-log(1i*x.*(1i*x+1))','x','C','p');
end
if all(p==0)
 pf=inline('-log(1i*x.*sign(2*C-1))','x','C','p');
end
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
%ak=goldsearch_spec(Kf,pf,anint,1e-6,k,Tau,r(:)',s0(:)',param,C,p);
%mean(ak)
ak=9.71*Ik;
z=x*Ik-1i*Ix*ak;
P=exp(-r.*Tau).*(we'*real(exp(Kf(z-1i*Ix*p,Ix*Tau,Ix*r,Ix*s0,param)+pf(z,Ix*C,Ix*p)+(-1i*z.*(Ix*k)))))/pi;
%P(~CC)=P(~CC)-(exp(s0(~CC))-exp(k(~CC)-r(~CC).*Tau(~CC)));

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

function [K,anint,anint2]=HestonKf(z,T,r,s0,param)
%function [K,anint]=HestonKf(z,T,r,s0,param)
%
% log characteristic function for the Heston model
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
%  v0=param(1);    : Initial vol.
%  kap=param(2);   : Mean reversion
%  th=param(3);    : Mean level
%  s=param(4);     : Vol diffusion.
%  rho=param(5);   : Correlation
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006 Magnus Wiktorsson
 
v0=param(1);    % Initial vol.
kap=param(2);   % Mean reversion
th=param(3);    % Mean level
s=param(4);     % Vol diffusion.
rho=param(5);   % Correlation

z2=-1i;
%a=((1i*z)-(1i*z).^2)/2;
%b=(kap-rho*s*1i*z);
%c=s^2/2;
%g=sqrt(4*a*c+b.^2);
%D=-2*a.*(1-exp(-g.*T))./((g-b)+(g+b).*exp(-g.*T));
anint2=sort(roots([1 -(s^2-2*rho*s*kap)/(1-rho^2)/s^2 -kap^2/(1-rho^2)/s^2]))';
if nargout>=2
 y=max(T(:));
 sigma=s;
 kappa=kap;
 a2=roots([(-rho^3*sigma^3*y^3+rho*sigma^3*y^3) (3*kappa*rho^2*sigma^2*y^3-rho*sigma^3*y^3-kappa*sigma^2*y^3+6*rho^2*sigma^2*y^2-6*sigma^2*y^2) (-3*kappa^2*rho*sigma*y^3+kappa*sigma^2*y^3-12*kappa*rho*sigma*y^2+6*sigma^2*y^2-24*rho*sigma*y) kappa^3*y^3+6*kappa^2*y^2+24*kappa*y+48]);
 anint=[max(a2(a2<=0)) min(a2(a2>0))];
end 
d=sqrt((kap-rho*s*1i*z).^2+s^2*(1i*z-(1i*z).^2));
CC=kap*th/s^2*((kap-rho*s*1i*z).*T-2*log(((kap-1i*rho*s*z)+d.*coth(d.*T/2))./(d./sinh(d.*T/2))));
%D=(1-exp(-d.*T)).*((1i*z).^2-1i*z)./((kap-1i*rho*s*z).*(1-exp(-d.*T))+d.*(exp(-d.*T)+1));
DD=((1i*z).^2-1i*z)./((kap-1i*rho*s*z)+d.*coth(d.*T/2));
K=1i*z.*(s0+r.*T)+CC+DD*v0;



