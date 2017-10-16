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

function [U,XXP,TT]=BSamPutUII_FGL(S,K,T,r,sigma)
% function [U,XXP,TTP]=BSamPutUII_FGL(S,K,T,r,sigma);
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
%   U: Black Scholes American put option price
%   XXP: optimal excise level
%   TTP: timepoints for XXP
%
% (C) 2014 Magnus Wiktorsson
[U,XXP,TT]=fftAP(S,@BSKf,K,T,r,sigma);
function [U,XXP,TT]=fftAP(SS,Kf,K0,T,r,par)
n=100;
TT=linspace(0,1,n+1);
TT=T*TT;
load glax1000
alpha=0.010422;
z=alpha+1i*x;
CC=exp((z+1)*log(K0))./z./(1+z);
dt=TT(2)-TT(1);
G=exp(Kf(1i*z,dt,r,0,par)-r*dt);
KKfr=Kf(1i*z,dt,r,0,par)-r*dt;
XXP=K0*ones(1,n+1);
for k=n:-1:1
CC=G.*CC;
ylast=log(XXP(n+1-k));
f=@(y) we'*real(exp(-y*z).*(CC+dt*K0*r*(exp(z*ylast+KKfr)-exp(z*y))./z./(KKfr+z*(ylast-y))))/pi-max(K0-exp(y),0);
[y]=binsearch(f,0.95*ylast,ylast,1e-12,k,n);
XXP(n+1-k+1)=exp(y);
CC=CC+dt*K0*r*(exp(z*ylast+KKfr)-exp(z*y))./z./(KKfr+z*(ylast-y));
end
IC=SS>XXP(end);
IS=ones(size(SS(IC)));
U=K0-SS;
U(IC)=real(we'*(exp(-z*log(SS(IC))).*(CC*IS)))/pi;




function [x]=binsearch(f,a,b,e,k,n)
fa=f(a);
fb=f(b);
N=max(ceil((log(b-a)-log(e))/log(2)));
if any((fa.*fb)>0)&(k>(0.9*n))
    k
   error('Same sign of f at endpoints') % will never happen (I hope)
else   
 for i=1:N
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
