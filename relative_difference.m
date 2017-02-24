% reldiff = relative_difference(left, right)
%
% compute the relative difference between two sets of numbers.
%
% copyright 2015, 2016 Daniel Brake
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

function reldiff = relative_difference(left, right)


reldiff = zeros(size(left));

indices_equal = (left-right) == 0;
indices_unequal = ~indices_equal;



reldiff(indices_unequal) = abs(left(indices_unequal) - right(indices_unequal)) ./ max(abs(left(indices_unequal)),abs(right(indices_unequal)));

end
