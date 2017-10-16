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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.%
% Note that we will use the European solver before and after the dividend
% Note also that we work with function values, and therefore a simple max
% is sufficient at the time of the dividend
%
function U=AmericanDDRBF1D(S,K,T,r,sig,payoff,phi,ep,xc,M,op,varep,alpha,D)
% Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
% Assuming a vector of proportional dividend times alpha and values D
    
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
  % Divide up the time steps between the intervals
  %
  M1 = round((1-alpha)*M);
  M2 = round(alpha*M);
  %
  % THE FOLLOWING PART IS FOR THE TIME FROM MATURITY BACKWARDS TO THE DIVIDEND
  %
  %---------------------------------------------------------------------------
  [k,beta0,beta1,beta2]=BDF2coeffs((1-alpha)*T,M1);
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
  % Initial condition. Note that boundaries are at the next time
  % 
  u0 = payoff((1-D)*xc,K,r,0);
  
  %figure(1),clf
  %plot(xc,u0,'k-'), hold on
  
  rhs(:,1) = u0;
  tn = k(1);
  rhs(end,1) = payoff((1-D)*x_max,K,r,tn);
  %
  % Time-stepping loop
  %
  for n = 1:M1
    %
    % New solution value from linear system
    %
    u = Uc\(Lc\rhs);
    %
    % Prepare the right hand side for the next step
    %
    nextstep = min(M1,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    rhs(1) = 0;
    rhs(end) = payoff((1-D)*x_max,K,r,tn);
    u0 = u;
  end
  
  %plot(xc,u,'b-')
  
  %--------------------------------------------------
  % Prepare for the next PDE problem
  [k,beta0,beta1,beta2]=BDF2coeffs(alpha*T,M2);
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
  % Now we need to check which ones fall below the payoff in x instead of y
  %
  u0 = max(u,payoff(xc,K,r,0));
  
  %plot(xc,u0,'r--')
  
  rhs(:,1) = u0;
  tn = k(1);
  rhs(end,1) = payoff(x_max,K,r,tn); % This is assuming that we had an
                                     % American issue
  %
  % Time-stepping loop
  %
  for n = 1:M2
    %
    % New solution value from linear system
    %
    u = Uc\(Lc\rhs);
    %
    % Prepare the right hand side for the next step
    %
    nextstep = min(M2,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead, Restarting time
                             % from zero        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    rhs(1) = 0;
    rhs(end) = payoff(x_max,K,r,tn);
    u0 = u;
  end
  %plot(xc,u,'b.-')
  %
  % Compute the solution or some derivative at the evaluation points
  %
  U = Ae*u;
  %
  % Compare with the reference values
  %
  U = U*fac;

