% function this = dehomogenize(this)
%
% dehomogenize a tropical curve.
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



function this = dehomogenize(this)

if ~(this.is_homogenized)
	error('cannot dehomogenize a non-homogenized curve')
end
num_vars = length(this.variables);

hom_index = 0;
for ii = 1:length(this.variables)
	if strcmp(this.variables{ii},this.homogenizing_variable)
		hom_index = ii;
	end
end
if hom_index==0
	error('unable to find homogenizing variable ''%s'' in set of variables',this.homogenizing_variable);
end


tmp_rays = zeros(length(this.variables)-1, size(this.rays,2));
ind = 1:num_vars;
for ii = 1:size(this.rays,2)
	tmp_rays(:,ii) = this.rays(ind~=hom_index,ii) - this.rays(hom_index,ii);
end





tmp_pts = zeros(length(this.variables)-1, size(this.aux_data.intersection_points,2));
ind = 1:num_vars;
for ii = 1:size(this.aux_data.intersection_points,2)
	h = this.aux_data.intersection_points(hom_index,ii);
	if abs(h)<1e-8
		h = 0;
	end
	
	tmp_pts(:,ii) = this.aux_data.intersection_points(ind~=hom_index,ii) ./ h;
end




tmp_vars = {};
for ii = 1:length(this.variables)
	if ii~=hom_index
		tmp_vars{end+1} = this.variables{ii};
	end
end

for ii = 1:length(this.variables)
	sym(this.variables{ii});
end


for ii = 1:size(this.curve_system,1)
	this.curve_system{ii,2} = subs(this.curve_system{ii,2}, genvarname(this.homogenizing_variable),1);
end



this.is_homogenized = false;
this.rays = tmp_rays;
this.aux_data.intersection_points = tmp_pts;
this.variables = tmp_vars;



for ii = 1:length(this.aux_data.paths)
	this.aux_data.paths(ii).dehomogenize(hom_index);
end





end
