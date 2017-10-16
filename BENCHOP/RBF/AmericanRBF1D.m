 %Copyright (C) 2015 Elisabeth Larsson

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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.% -dt*Q in the right hand side, but what about coeffs? Written in a form
% such that L only has dt coeff.
% Q = (e*C-e*(d1*xc(:,1)+d2*xc(:,2))/2)./(u(:,n-1)+e-q); Föregående
% tidsnivå används.
% q = (xc(:,1)+xc(:,2))/2-K; (No maximum taken)
% d1,d2 är dividends. Here we have put instead of call, so different term
% and different bcs.
% Asymptotic for put?
% C = r*K;
function U=AmericanRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,epen)
% Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
  
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
  K0 = K;
  K = 1;
  %
  % Evaluation grid points (uniform)
  %
  xe = S/K0;
  xc = xc/K0;
  N = length(xc);
  x_max = xc(end);
  p = str2num(op);
  fac = K0^(1-p);

  x_min = 0;
  sigma2 = sig^2;

  [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
  %
  % If variable epsilon is requested
  %
  if strcmp(lower(varep),'yes')
    d(1)=xc(2)-xc(1);
    d(2:N-1,1) = 0.5*(xc(3:N)-xc(1:N-2));
    d(N) = xc(N)-xc(N-1);
    d = d/max(d);
    % We let the scaling be between ep/Q and ep
    % ep(dmax)=ep/Q, ep/Q/(dmin+alpha)=ep, dmin+alpha = 1/Q
    Q=16;
    ep0=ep;
    alpha = max(0,1/Q-min(d));
    ep = ep0/Q./(alpha+d);
  end  
  %
  % Diagonal matrices with coordinates for coefficients
  %
  X = spdiags(xc(:,1),0,N,N);
  %
  % Distance matrix to build interpolation and collocation matrices
  %
  rc = xcdist(xc,xc,1);
  re = xcdist(xe,xc,p);
  %
  % Interpolation and operator matrices
  %
  A0 = RBFmat(phi,ep,rc,'0');
  Ax = RBFmat(phi,ep,rc,'1',1);
  Axx = RBFmat(phi,ep,rc,'2',1);
  Ae = RBFmat(phi,ep,re,op);
  %
  % Build the Black-Scholes operator matrix
  %
  L = (sigma2/2*X.^2*Axx + r*X*Ax  - r*A0)/A0;
  %
  clear Ax Axx 
  %
  % The differentation matrix A*A^-1
  % 
  Ae = Ae/A0;
  %
  % Time-stepping matrix I-beta_0*L in the interior
  %
  C = eye(N) - beta0*L;
  %
  % Unit one in the rows corresponding to boundary conditions
  %
  C(1,:) = 0;
  C(end,:) = 0;
  C(1,1) = 1;
  C(end,end) = 1;
  [Lc,Uc] = lu(C);
  %
  % Prepare for the non-linear penalty term
  %
  %epen = 8e-7;
  q = K-xc; %(No maximum taken)
  %
  % Initial condition. Note that boundaries are at the next time
  % 
  u0 = payoff(xc,K,r,0); 
  P = (epen*r)./(u0+epen-q); % Divided by K above and below;
  rhs(:,1) = u0 + beta0*P; 
  tn = k(1);
  rhs(1) = payoff(x_min,K,r,tn);
  rhs(end,1) = payoff(x_max,K,r,tn);
  %
  % Time-stepping loop
  %
  for n = 1:M
    %
    % New solution value from linear system
    %
    u = Uc\(Lc\rhs);
    %
    % Prepare the right hand side for the next step
    %
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    P = (epen*r)./(u+epen-q); % Divided by K above and below
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0 + beta0*P;
    % boundary conditions
    rhs(1) = payoff(x_min,K,r,tn);
    rhs(end) = payoff(x_max,K,r,tn);
    u0 = u;
  end
  %
  % Compute the solution or some derivative at the evaluation points
  %
  U = Ae*u;
  %
  % Compare with the reference values
  %
  U = U*fac;

