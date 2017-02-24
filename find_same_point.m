%  indices = find_same_point(point_collection, test_point, tolerance, how_many)
%
%  finds indices of test_point in point_collection, columnwise, up to
%  tolerance.  for coordinates which are large, uses relative difference of
%  the coordinates instead of the absolute difference.
%
% if how_many is specified, this call is passed on to Matlab's find()
% command, to only find that many (if they exist).
%
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

function indices = find_same_point(point_collection, test_point, tolerance, how_many)



R = repmat(test_point,[1 size(point_collection,2)]);


bigind = or(abs(point_collection)>1,abs(R)>1);

testme = abs(point_collection - R);

testme(bigind) = relative_difference(point_collection(bigind),R(bigind));

% ind represents the indices of the intersection point which the lower
% point tracked TO.
if nargin<=3
	indices = find(sum(testme > tolerance,1)==0);
else
	indices = find(sum(testme > tolerance,1)==0, how_many);
end

end
