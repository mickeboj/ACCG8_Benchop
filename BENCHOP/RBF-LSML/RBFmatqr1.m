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

function [A,Ax,Axx]=RBFmatqr1(phi,ep,xe,xc,p)
%
% Should select if to use standard method or to use qr  
% However, we need reuse. Nargout determines the number of derivatives.
%  
  %R = (xc(end)-xc(1))/2; % Scaling compared with ref domain [-1 1]
  %if (ep*R >= 2) % Use direct method
  if (nargin==4 & nargout > 1 | nargin==5 & p>0)
    r = xcdist(xe,xc,1); 
    dim = 1;
    if (nargin==4)
      A = RBFmat(phi,ep,r,'0');
      if (nargout>1)
        Ax = RBFmat(phi,ep,r,'1',dim);
      end  
      if (nargout >2)
        Axx = RBFmat(phi,ep,r,'2',dim);
      end
    else
      op = num2str(p);
      A = RBFmat(phi,ep,r,op,dim);
    end  
  else
    r = xcdist(xe,xc,0);
    A = RBFmat(phi,ep,r,'0');
  end  







