% [intersection_points] = compute_intersection_points(b_input, aux_data, reality_string)
%
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

function [intersection_points] = compute_intersection_points(b_input, aux_data, reality_string)

if exist('previous_intersection_points.mat','file')
	load('previous_intersection_points.mat');
	if issame_system(previous_b_input,b_input)
		display('reusing previous intersection points');
		intersection_points = previous_intersection_points;
		return;
	end
end



	

system_variables = b_input.variable_group;
num_vars = numel(system_variables);


% for each variable, compute the intersection points with the corresponding axis.
intersection_points = [];
for kk = 1:num_vars
	display(sprintf('slicing %s=0',system_variables{kk}));
	local_intersection_points = coordinate_slice(b_input, aux_data.initial_start_points, aux_data.random_line_coefficients, 0, kk, reality_string);
	intersection_points = [intersection_points local_intersection_points];
end



saveme.previous_intersection_points = intersection_points;
saveme.previous_b_input = b_input; %#ok 
save('previous_intersection_points.mat','-struct','saveme')
clear saveme

end
