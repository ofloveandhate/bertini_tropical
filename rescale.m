% rescale a single point, or a matrix containing points in a
% direction.
%
% if you wish to rescale a point, just pass it in, and it will rescale based on the largest coordinate.
%
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

function result = rescale(point,dimension_for_rescale)




if nargin < 2
	dimension_for_rescale = 1;
end



if ndims(point)==2
	if sum(size(point)==1)==1
		[~,i] = max(abs(point));
		result = point/point(i);
		return
	end
	
	if dimension_for_rescale==0
		error('must specify a dimension for rescaling');
	end
	if dimension_for_rescale==1
		[m,i] = max(abs(point),[],dimension_for_rescale);
		ind = sub2ind(size(point),i,1:size(point,2));
		result = point ./ abs(repmat(point(ind),[size(point,1) 1]));
	elseif dimension_for_rescale==2
		[m,i] = max(abs(point),[],dimension_for_rescale);
		ind = sub2ind(size(point),1:size(point,1),i');
		result = point ./ abs(repmat(point(ind)',[1 size(point,2)]));
	else
		error('specified dimension exceeds dimension of object to rescale');
	end
else
	error('rescaling for an object of %i dimensions is not implemented',ndims(point));
end


end
