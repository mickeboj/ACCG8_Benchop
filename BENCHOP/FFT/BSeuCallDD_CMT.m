 %Copyright (C) 2015 Erik Lindström

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

function [U,errRel]=BSeuCallDD_CMT(S,K,T,r,sigma,D,alpha)
S=S(:)';
U=fftEU(@BSDDKf,K,S,r,T,[sigma D alpha*T]);


function [P0]=fftEU(Kf,K0,S0,r,T,par)

sig=par(1);
D=par(2); 
tau=par(3);

% tic, 
m=2^8;
Zmax=2*pi/sqrt(sig^2*T)/40*m;
dz=Zmax/m;
z=(0:(m-1))*dz;
dx=2*pi/dz/m;
b=m/2*dx;
x=b-(0:(m-1))*dx;
S=exp(x+log(K0));
alpha=30;
w=ones(1,m);w([1 m])=1/2;
CC=exp(1i*b*z)./(1i*z+alpha)./(1i*z+alpha+1).*exp(Kf(z-1i*(alpha+1),T,r,0,par)-r*T);
P=dz*exp(alpha*x+x+log(K0)).*real(fft(CC.*w))/pi;
F=griddedInterpolant(fliplr(S),fliplr(P),'spline');
P0=F(S0);
% toc


% % Compute relative error
% P0_BSDD=BSDDCall_exact(S0,K0,T,r,par);
% errRel=(P0-P0_BSDD)./P0_BSDD

% % Plot solution
% SGrid=60:1:160;
% plot(SGrid,F(SGrid),SGrid,BSDDCall_exact(SGrid,K0,T,r,par),'--');
% legend('Carr-Madan FFT','Exact')
% xlabel('S')
% ylabel('Price')


function [K,anint]=BSDDKf(u,T,r,s0,par)
% log characteristic function for the Black-Scholes model
% discrete divedend at tau
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
%     s: par(1) BM-vol
% delta: par(2) relative dividend
%   tau: par(3) time of dividend 
%
% ---------------output-------------------
% 
%      K : log characteristic function 
%   anint: strip of analyticity for current parameter values 

% (C) 2006,2014 Magnus Wiktorsson

s=par(1); % BM-vol
delta=par(2); % relative dividend
tau=par(3); % time of dividend 
K=1i*u.*(s0+T.*(r-s.^2/2))-T.*u.^2/2*s.^2+1i*u.*(tau<=T)*log(1-delta);
anint=[-inf,inf];




function [CallBSDD]=BSDDCall_exact(S,K,T,r,par)
% function [CallBSDD]=BSCall_exact(S,K,T,r,par);
%
% ----------input---------------------------
%
%       S: intial/current stock price
%       K: strike price
%       T: time to maturity
%       r: continously compounded short rate
%
%     par: parameters lack Scholes discrete dividend model
%
%     sigma: par(1) BM-vol
%         D: par(2) relative dividend
%       tau: par(3) time to dividend 
%
%   
%  Note that for S,K,T and r that they should either be 1-dim
%  or have matching dimensions
%
%  If tau>T then the price collapses to the standard Black Scholes price
%
% ----------output---------------------------
%
%   CallBS: Black Scholes discrete dividend European call option price
%
% (C) 2014 Magnus Wiktorsson

     sigma=par(1); % BM-vol
         D=par(2); % relative dividend
       tau=par(3); % time to dividend 


CallBSDD=S.*(1-D.*(tau<T)).*normcdf((log(S./K)+(r+sigma.^2/2).*T+log(1-D).*(tau<=T))./(sigma.*sqrt(T)))-K.*exp(-r.*T).*normcdf((log(S./K)+(r-sigma.^2/2).*T+log(1-D).*(tau<=T))./(sigma.*sqrt(T)));
