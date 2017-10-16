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

function U = EuropeanVegaLSML(S,K,T,r,sig,payoff,phi,epf,epc,xf,xc,xls,M)
%
% Main program for all the different variants of the method 
%
% phi is the RBF 'mq', 'gs',...
% ep(1:L) is the (constant) shape parameter for each grid
%         Note that for a one grid method ep is just given as a scalar 
% M is the number of time steps
% N(1:L) is the number of center points for each grid
% Ns is the number of uniform points where the solution is evaluated
% rg(1:2) is the range over which we are computing, typiclly [0 4]
% ctype is 'cheb' or 'uni' for the distribution of node points
% etype has the same function for the points where equations are enforced
% method{1:L} is 'qr' or 'dir' for RBF-QR or RBF-Direct for each grid
%             note only the lsml version allows different
% show is 'yes' or 'no' and tells if plots of the errors should be made
% col is the color to use for the error plots.
% Nls is _optional_. If present it indicates that least squares should be
% used.
%-----------------------------------------
% The parameters are given unscaled, but in the solver we always scale
% down x and U so that K=1.
  op='0';
  
  K0 = K;
  K = 1;
  %
  % Scale the nodes.
  %
  Nc = length(xc);
  xc = xc/K0;
  Nb = 2;
  xb = [xc(1);xc(end)];
  Nf = length(xf);
  xf = xf/K0;
  Nls = length(xls);
  xls = xls/K0;
  xe = S/K0;
  %
  p = str2num(op);
  fac = K0^(1-p);
  %  
  % BDF2-coefficients to have constant matrices. 
  %  
  [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
  %
  % Initiate the matrices that we need for time-stepping (both levels)
  %
  X = spdiags(xls(:,1),0,Nls,Nls);
  %
  if (Nf >0)
    A0 = RBFmatqr1(phi,epf,xf,xf);
    [Af,Axf,Axxf] = RBFmatqr1(phi,epf,xls,xf);
    Lf = (sig^2/2*X.^2*Axxf + r*X*Axf  - r*Af)/A0;
    If = Af/A0;
    Cf = If-beta0*Lf;
    Ef = RBFmatqr1(phi,epf,xe,xf,p)/A0;
    Axxf = Axxf/A0; 

    Ifb = If(:,[1 Nf]); If = If(:,2:Nf-1);
    [Qf,Rf]=qr(Cf(:,2:Nf-1));
    Cf = Cf(:,[1 Nf]);
  end

  %
  A0 = RBFmatqr1(phi,epc,xc,xc);
  [Ac,Axc,Axxc] = RBFmatqr1(phi,epc,xls,xc);
  Lc = (sig^2/2*X.^2*Axxc + r*X*Axc  - r*Ac)/A0;
  Ic = Ac/A0;
  Cc = Ic-beta0*Lc;
  Icb = Ic(:,[1 Nc]); Ic = Ic(:,2:Nc-1);
  [Qc,Rc]=qr(Cc(:,2:Nc-1));
  Cc = Cc(:,[1 Nc]);  
  Ec = RBFmatqr1(phi,epc,xe,xc,p)/A0;
  Axxc = Axxc/A0;
  %
  % Evaluation matrices NOTE:Bug in the sense that I may be requesting
  % only Ax here!!
  %
  %
  % Initial condition
  %
  u0 = payoff(xls,K,r,0);
  phi0 = zeros(size(u0));
  %
  % Initial rhs. Subtraction of the boundary part comes in each step
  %
  rhs0 = u0;  
  phirhs = phi0;
  %
  % The time-stepping loop
  %
  res = 0; timeres=0;
  for step = 1:M
    tn = sum(k(1:step));
    gn = payoff(xb,K,r,tn); % For phi, gn=0;
    rhs_c = rhs0 - Cc*gn;
    %
    % Solve for interior nodes at coarse level
    %
    bc = Qc'*rhs_c;
    uc = Rc\bc;
    %
    % Compute residual (rhs for fine level). Note that g=0 for this step.
    %
    % Rc*uc=bc Q*bc= rhs_c. It does make sense in both forms, but gives 
    % different results.
    %rhs_f = rhs_c - Qc*bc;
    if (Nf > 0)
      rhs_f = rhs_c - Qc*Rc*uc;
      %diff=max(abs(rhs_f0-rhs_f))
      %
      % Solve for interior nodes at fine level
      %
      bf = Qf'*rhs_f;
      uf = Rf\bf;
      %check=max(abs(uf))
      % uf = 0*uf;
      %
      % Evaluate two-level solution at least squares points
      %
      u = If*uf+ Ic*uc + Icb*gn;
    else
      u = Ic*uc + Icb*gn; 
    end
    %
    % Compute right hand side for next step
    %
    nextstep = min(step+1,M);
    rhs0 =  beta1(nextstep)*u - beta2(nextstep)*u0;
    %
    % Move the values for the next step
    %
    u0 = u;
    %
    % The computation of phi is one step after. We use the newly computed u
    % Think about the times. Perhaps write it down if it seems to work.
    % 
    % First evaluate the second derivative for the right hand side
    uxx = Axxc*[gn(1);uc;gn(end)];
    if (Nf>0)
      uxx = uxx + Axxf*[0;uf;0];
    end
    phirhs = phirhs + beta0*sig*K0*X.^2*uxx;
    %
    % Solve for phi
    %
    bc = Qc'*phirhs;
    phic = Rc\bc;
    phi = Ic*phic;
    %
    % Compute residual (rhs for fine level). Note that g=0 for this step.
    %
    % Rc*uc=bc Q*bc= rhs_c. It does make sense in both forms, but gives 
    % different results.
    %rhs_f = rhs_c - Qc*bc;
    if (Nf > 0)
      phirhs_f = phirhs - Qc*Rc*phic;
      %diff=max(abs(rhs_f0-rhs_f))
      %
      % Solve for interior nodes at fine level
      %
      bf = Qf'*phirhs_f;
      phif = Rf\bf;
      %check=max(abs(uf))
      % uf = 0*uf;
      %
      % Evaluate two-level solution at least squares points
      %
      phi = phi + If*uf;
    end
    %
    % Prepare next rhs for phi
    %
    phirhs = beta1(nextstep)*phi-beta2(nextstep)*phi0;
    %
    % Move the values for next step
    %
    phi0 = phi;
  end  
  %
  % Compute the values at the final time
  %
  U = Ec(:,2:end-1)*phic;
  if (Nf>0)
    U = U + Ef(:,2:end-1)*uf;
  end

