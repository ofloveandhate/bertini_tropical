% this = rescale(this)
%  Rescale a tropical curve.
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



function this = rescale(this)


tmp_pts = zeros(length(this.variables), size(this.aux_data.intersection_points,2));
for ii = 1:size(this.aux_data.intersection_points,2)
	tmp_pts(:,ii) = rescale(this.aux_data.intersection_points(:,ii));
end


this.aux_data.intersection_points = tmp_pts;


for ii = 1:length(this.aux_data.paths)
	this.aux_data.paths(ii).rescale();
end

for ii = 1:length(this.aux_data.corresponding_centers)
	this.aux_data.corresponding_centers{ii} = rescale(this.aux_data.corresponding_centers{ii});
end



end
