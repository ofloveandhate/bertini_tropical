% [coefficient] = cauchy_integral(path,k)
%
% path is a path around the origin, and k is the power of the denomenator
% given a path around the origin (path), integrate it using the trapezoidal
% rule.
% 
% copyright 2016 Daniel Brake
% University of Notre Dame
% Applied and Computational Mathematics and Statistics
% danielthebrake@gmail.com
%
%  Bertini (TM) is a registered trademark.
%
%     This file is part of Bertini_tropical.
%
%     Bertini_tropical is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     Bertini_tropical is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function coefficient = cauchy_integral(path, radius, cycle_number, k)

nsegments = length(path)-1;
theta = (radius^(1/cycle_number)*exp( 1i * linspace(0,2*pi,nsegments+1) ) ).^k;
coefficient = 1/(nsegments) * sum(path(1:end-1)./theta(1:end-1));


end


