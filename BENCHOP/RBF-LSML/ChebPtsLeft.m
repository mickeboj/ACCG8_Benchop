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

function x=ChebPtsLeft(rg,N,ell);
  % Normal right cheb interval has length L=1 and starts at a=-1
  L = -cos(ell*pi); 
  a = -L ;  
  theta=(ell*pi:-(ell-1/2)*pi/(N-1):pi/2)';
  x = rg(1) + (rg(2)-rg(1))/L*(cos(theta)-a); 

