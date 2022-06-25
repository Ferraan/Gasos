## Copyright (C) 2022 Ferran
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} cp (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Ferran <ferran@SuperPC>
## Created: 2022-06-25

function cp = cp (T0, Tf,a0,a1,a2,a3,a4)
  fun = @(T) a0+a1*T+a2*T.^2+a3*T.^3+a4*T.^4;
  cp=integral(fun,T0,Tf)/(Tf-T0);
  
endfunction
