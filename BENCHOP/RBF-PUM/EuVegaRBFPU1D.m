 %Copyright (C) 2015 Victor Shcherbakov

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

function [Phi]=EuVegaRBFPU1D(S,K,T,r,sig,payoff,phi,ep,xc,M,np,op)
% Compute the value of Vega of a European option using uniform RBF-PUM approximation
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% N is the number of points in one spatial dimension
% M is the number of time steps
% Ne is the number of evaluation points in one dimensions
% x_max is the size of the computational domain in each dimension
% xc is computational points
% S evaluation points, points we want to know the price at
% K strike
% T maturity time
% r risk free interest rate
% sig volatility
% payoff payoff function
% op number of derivative
% np number of partitions
%
%  (C) Elisabeth Larsson & Victor Shcherbakov


% Define domain parameters
xeval = S;
Ne = length(S);

N = length(xc);
x_min = xc(1);
x_max = xc(end);

h = x_max/(N-1);

% Define weights on BDF2
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);

%
% Compute the size of the partitions
%
H = (x_max-x_min)/np;
x = H/2:H:x_max-H/2;
overlap = 3*h;
R = H/2 + overlap;
%
% Compute J, the maximum number of overlapping partitions at one point
%
n = floor(R/(H/2));
if rem(R,H/2)>0
  J = n+1;
else
  J = n;
end  

if (J~=2)
  error(['The chosen grid size and partition size leads to more than two ' ...
         'partitions overlapping or none']);
end  
%
% Check which partitions each node falls into
%
NodePart = zeros(N,J) - 1;
NodePart(:,1) = floor(xc/H) + 1;
NodePart(end,1)=np; %(Last point gets exactly np+1)

pos = find (rem(xc,H)<=overlap & NodePart(:,1) >1);
NodePart(pos,2) = NodePart(pos,1)-1;
pos = find (rem(xc,H)>=H-overlap & NodePart(:,1) <np);
NodePart(pos,2) = NodePart(pos,1)+1;
NodePart(end,2) = -1; % To not mistakenly place it in one more 
%
% Check which nodes fall into which partition
%
for i=1:np
  indlocpts{i} = find(NodePart(:,1)==i | NodePart(:,2)==i);
  Ni(i)=length(indlocpts{i});
end  
%
% Map global indices to local indices
%
LocPos = spalloc(N,N,sum(Ni));
for i=1:np
  LocPos(indlocpts{i},i) = 1:Ni(i);
end  

%
% Compute the weight functions for each partition using Shepard's method
%
for i = 1:np
  s{i}=zeros(Ni(i),1);
  sx{i}=zeros(Ni(i),1);
  sxx{i}=zeros(Ni(i),1);
end  
for i=1:np
  [psi{i},psix{i},psixx{i}] = RBFPUweight1D(xc(indlocpts{i}),x(i),R);
  s{i} = s{i} + psi{i};
  sx{i}= sx{i} + psix{i};
  sxx{i} = sxx{i} + psixx{i};
  if (i>1)
    posi = find(any(NodePart(indlocpts{i},:)==i-1,2));
    posm = LocPos(indlocpts{i}(posi),i-1);
    s{i-1}(posm) = s{i-1}(posm) + psi{i}(posi); 
    sx{i-1}(posm) = sx{i-1}(posm) + psix{i}(posi); 
    sxx{i-1}(posm) = sxx{i-1}(posm) + psixx{i}(posi); 
  end
  if (i<np)
    posi = find(any(NodePart(indlocpts{i},:)==i+1,2));
    posm = LocPos(indlocpts{i}(posi),i+1);
    s{i+1}(posm) = s{i+1}(posm) + psi{i}(posi); 
    sx{i+1}(posm) = sx{i+1}(posm) + psix{i}(posi); 
    sxx{i+1}(posm) = sxx{i+1}(posm) + psixx{i}(posi); 
  end  
end  


% Define partition of unity

sw=zeros(N,1);
for i=1:np
  s2 = s{i}.*s{i};
  s3 = s2.*s{i};
  W{i} = psi{i}./s{i};
  Wx{i} = psix{i}./s{i} - psi{i}.*sx{i}./s2;
  Wxx{i} = -2*psix{i}.*sx{i}./s2 + psixx{i}./s{i} + ...
           psi{i}.*(2*sx{i}.^2./s3 - sxx{i}./s2);
  sw(indlocpts{i}) =   sw(indlocpts{i}) + W{i};
end


%
% Define partition of unity for evaluation points
%
% Check which partitions each node falls into
%
NodePartE = zeros(Ne,J) - 1;
NodePartE(:,1) = floor(xeval/H) + 1;
pose = find(NodePartE(:,1)==np+1);
NodePartE(pose,1)=np; %(Last point gets exactly np+1)

pos = find (rem(xeval,H)<=overlap & NodePartE(:,1) >1);
NodePartE(pos,2) = NodePartE(pos,1)-1;
pos = find (rem(xeval,H)>=H-overlap & NodePartE(:,1) <np);
NodePartE(pos,2) = NodePartE(pos,1)+1;
NodePartE(pose,2) = -1; % To not mistakenly place it in one more 
%
% Check which nodes fall into which partition
%
for i=1:np
  Elocpts{i} = find(NodePartE(:,1)==i | NodePartE(:,2)==i);
  Nei(i)=length(Elocpts{i});
end  
eparts = find(Nei>0);
%
% Map global indices to local indices
%
ELocPos = spalloc(N,N,sum(Nei));
for i=eparts
  ELocPos(Elocpts{i},i) = 1:Nei(i);
end  
%
% Compute the weight functions for each partition using Shepard's method
%
for i = eparts
  s{i}=zeros(Nei(i),1);
end  
for i=eparts
  [psi{i},psix{i},psixx{i}] = RBFPUweight1D(xeval(Elocpts{i}),x(i),R);
  s{i} = s{i} + psi{i};
  if (i>1)
    posi = find(any(NodePartE(Elocpts{i},:)==i-1,2));
    if (length(posi)>0)
      posm = ELocPos(Elocpts{i}(posi),i-1);
      s{i-1}(posm) = s{i-1}(posm) + psi{i}(posi); 
    end
  end
  if (i<np)
    posi = find(any(NodePartE(Elocpts{i},:)==i+1,2));
    if (length(posi)>0)
      posm = ELocPos(Elocpts{i}(posi),i+1);
      s{i+1}(posm) = s{i+1}(posm) + psi{i}(posi); 
    end
  end  
end  


% Define partition of unity
sw=zeros(Ne,1);
for i=eparts
  We{i} = psi{i}./s{i};
  sw(Elocpts{i}) =   sw(Elocpts{i}) + We{i};
end

% Form local RBF matrices
L = spalloc(N,N,sum(Ni.^2));
Vxx = spalloc(N,N,sum(Ni.^2));
for i = 1:np
    
    xloc = xc(indlocpts{i});
    rc = xcdist(xloc,xloc,1);
    A0  = RBFmat(phi,ep,rc,'0');
    Ax  = RBFmat(phi,ep,rc,'1');
    Axx = RBFmat(phi,ep,rc,'2');

    Wi = spdiags(W{i},0,Ni(i),Ni(i));
    Wix = spdiags(Wx{i},0,Ni(i),Ni(i));
    Wixx = spdiags(Wxx{i},0,Ni(i),Ni(i));
    
    X = spdiags(xloc,0,Ni(i),Ni(i));
    
    D0 = Wi*A0; 
    Dx = Wi*Ax + Wix*A0;
    Dxx = Wi*Axx + 2*Wix*Ax + Wixx*A0; 

    B = (r*X*Dx + 0.5*sig^2*X.^2*Dxx - r*D0)/A0;

    if (find(eparts==i))
      xe = xeval(Elocpts{i});
      re = xcdist(xe,xloc,op);
      Wei = spdiags(We{i},0,Nei(i),Nei(i));
      E{i}  = Wei*RBFmat(phi,ep,re,op)/A0;
    end

    L(indlocpts{i},indlocpts{i}) = L(indlocpts{i},indlocpts{i}) + B;
    
    % Matrix of 2nd derivative, needed to compute Vega
    Vxx(indlocpts{i},indlocpts{i}) = Vxx(indlocpts{i},indlocpts{i}) + X.^2*Dxx/A0;
    
end


% Time-stepping matrix I-beta_0*L in the interior
C = speye(N) - beta0*L;

% Unit one in the rows corresponding to boundary conditions
C(1,:) = 0;
C(end,:) = 0;
C(1,1) = 1;
C(end,end) = 1;

[Lc,Uc] = lu(C);

% Initial condition. Note that boundaries are at the next time
u0 = payoff(xc,K,r,0); 
phi0 = zeros(size(u0));
phirhs = phi0; 
rhs(:,1) = u0;
tn = k(1);
rhs(end,1) = payoff(x_max,K,r,tn);


% Time-stepping loop
for n = 1:M

    % New solution value from linear system
    u = Uc\(Lc\rhs);

    % Prepare the right hand side for the next step
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead        
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    % boundary conditions
    rhs(1) = 0;
    rhs(end) = payoff(x_max,K,r,tn);
    u0 = u;
    
    % The computation of phi is one step after. We use the newly computed u
    % Think about the times. Perhaps write it down if it seems to work.
    % 
    phirhs = phirhs + beta0*sig*(Vxx*u);
    phirhs(1)=0;
    phirhs(end)=0;
    phi = Uc\(Lc\phirhs);
    phirhs = beta1(nextstep)*phi-beta2(nextstep)*phi0;
    phi0 = phi;
end

Phi = zeros(Ne,1);

% Evaluate solution at needed (evaluation) points
for i=eparts
  Phi(Elocpts{i}) = Phi(Elocpts{i}) + E{i}*phi(indlocpts{i}); 
end 
 
end





