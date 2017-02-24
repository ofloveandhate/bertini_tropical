% dehomogenize a single point, or a matrix containing points in a
% direction.
%
% if you wish to dehomogenize a point, just pass it in, and it will assume
% the first coordinate is the homogenizing coordinate.
% if you wish to dehom a point, and the first coord is not the hom
% variable, specify its index as the third index
%
%
%
% copyright 2013 - 2016 Daniel Brake
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

function result = dehomogenize(point,dimension_for_dehom, index_of_dehom_coord)

if nargin < 3
	index_of_dehom_coord = 1;
end


if nargin < 2
	dimension_for_dehom = 1;
end


if length(point)<2
	error('point too short to dehomogenize');
end



if ndims(point)==2
	if sum(size(point)==1)==1
		result = point((1:end)~=index_of_dehom_coord)/point(index_of_dehom_coord);
		return
	end
	if nargin==1
		error('must specify a dimension for dehomogenization');
	end
	if dimension_for_dehom==1
		result = point( (1:end)~=index_of_dehom_coord,:) ./ repmat(point(index_of_dehom_coord,:),[size(point,1)-1 1]);
	elseif dimension_for_dehom==2
		result = point(:,(1:end)~=index_of_dehom_coord) ./ repmat(point(:,index_of_dehom_coord),[1 size(point,2)-1]);
	else
		error('specified dimension exceeds dimension of object to dehomogenize');
	end
else
	error('dehomogenization for an object of %i dimensions is not implemented',ndims(point));
end


end
