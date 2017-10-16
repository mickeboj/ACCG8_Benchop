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

function x=ChebPts(rg,N,ell);
% ell = 1 is the normal case
% smaller ell makes less stretching
% ell should not be smaller than or equal to 0.5  
  % Normal cheb interval has length L=2 and starts at a=-1
  L=2*cos((1-ell)*pi); 
  theta=(pi*ell:-pi*(2*ell-1)/(N-1):(1-ell)*pi)';
  sc = (rg(2)-rg(1))/L;
  x = rg(1) + sc*(cos(theta)-cos(theta(1))); % Make nodes in [0 D]


  
  % Normal clustering, -pi,0 If we want to change we do -pi*ell:(1-ell)*pi
  % Removing (1-ell)*pi in both ends
  
