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
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function U=EuropeanRBF2D(S,K,T,r,sig1,sig2,rho,payoff,phi,ep,xc,xtype,M,op,varep)
% Compute the value of a European option using uniform RBF approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
  
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
  x_max = max(max(xc));
  K0 = x_max;%x_max/4;
  %
  % Evaluation grid points
  %
  K = K/K0;
  xc = xc/K0;
  xe = S/K0;
  
  N = size(xc,1);
  %
  % The boundary is now 2-D and needs to be identified
  %
  b0 = find(xtype==1);
  bp = find(xtype==2);
  %bq = find(xtype==3);
  p = str2num(op);
  fac = K0^(1-p);

  x_min = 0;
  [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
  %
  % If variable epsilon is requested
  %
  if strcmp(lower(varep),'yes')
    %
    % For each point compute the distance to two of the nearest
    % neighbours. Use this to determine epsilon
    %
    d(1)=norm(xc(2,:)-xc(1,:));
    d(2:N-1,1) = min( sum((xc(2:N-1,:)-xc(1:N-2,:)).^2,2),...
                      sum((xc(3:N,:)-xc(2:N-1,:)).^2,2) );
    d(2:N-1,1) = sqrt(d(2:N-1,1));
    d(N) = norm(xc(N,:)-xc(N-1,:));
    d = d/max(d);
    % We let the scaling be between ep/Q and ep
    % ep(dmax)=ep/Q, ep/Q/(dmin+alpha)=ep, dmin+alpha = 1/Q
    Q=16;
    ep0=ep;
    alpha = max(0,1/Q-min(d));
    ep = ep0/Q./(alpha+d);
    epmin=min(ep)
    epmax=max(ep)
  end  
  %
  % Diagonal matrices with coordinates for coefficients
  %
  X = spdiags(xc(:,1),0,N,N);
  Y = spdiags(xc(:,2),0,N,N);
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
  Ay = RBFmat(phi,ep,rc,'1',2);
  Ayy = RBFmat(phi,ep,rc,'2',2);
  Axy = RBFmat(phi,ep,rc,'m2',[1 2]);
  Ae = RBFmat(phi,ep,re,op);
  %
  % 
  %
  [L0,U0]=lu(A0);
  %
  % Build the Black-Scholes operator matrix
  %
  I = speye(N);
  %
  L = ((1*(sig1^2/2*X.^2*Axx + rho*sig1*sig2*X.*Y*Axy + sig2^2/2*Y.^2*Ayy) + r*X*Ax ...
       + r*Y*Ay - r*A0)/U0)/L0;
  %
  clear Axx Ax Axy Ayy Ay 
  %
  % The differentation matrix A*A^-1
  % 
  Ae = (Ae/U0)/L0;
  %
  % Time-stepping matrix I-beta_0*L in the interior
  %
  C = I - beta0*L;
  %Ctest = cond(C)
  %
  % Unit one in the rows corresponding to Dirichlet boundary conditions
  %
  C(bp,:) = I(bp,:);

  C(b0,:) = I(b0,:);


  [Lc,Uc] = lu(C);
  %
  % Initial condition. Note that boundaries are at the next time
  %
  tn = k(1);
  u0 = payoff(xc(:,1),xc(:,2),K,r,0);

  rhs(:,1) = u0;
  rhs(bp,1) = payoff(xc(bp,1),xc(bp,2),K,r,tn);
  rhs(b0,1) = 0;
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
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    rhs(bp) = payoff(xc(bp,1),xc(bp,2),K,r,tn);
    rhs(b0) = 0;
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

