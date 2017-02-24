% slice a curve around a coordinate axis, slice it above and below, 
% and then walk around the t-origin, determining the individual paths
% converging to the t-origin.
%
% slice_value: the numeric value at which to slice.  slices at both plus
% and minus this value.  this should be determined previously by a critical
% point computation.
%
% local_isect: precomputed intersection points
%
% initial_start_points: precomputed witness points to move from.
%
% random_line_coefficients: the coefficients of the random linear slice
% (initial start points must lie on this line).  this is probably produced
% by a tracktype1 bertini run during NID.
%
% component_degree: the degree of the complex component.  determines an
% upper bound for the number of monodromy loops you might have to walk.
%
% rest: the definition of the system we are working with.
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

function [walked_paths,reject_connecting_points,slice_points, connections] = slice_and_monodromy_complex(slice_value, intersection_points, initial_start_points, random_line_coefficients, component_degree, b_input_copyme, variable_index, num_sample_points, unique_point_threshold)


b_input = copy(b_input_copyme);

system_variables = b_input.variable_group;
%switch directories to a temp directory.
dirname = sprintf('cauchy_%s',system_variables{variable_index});
if ~isdir(dirname)
	mkdir(dirname);
end
cd(dirname);


if ~isempty(dir('previous_paths.mat'))
	load('previous_paths.mat');
	if issame_system(previous_b_input,b_input) && (previous_slice_value==slice_value)
		display('reusing previous paths');
		walked_paths = previous_paths;
		reject_connecting_points = previous_reject_connecting_points;
		slice_points = previous_slice_points;
		connections = previous_connections;
		cd ..
		return;
	end
end 




% compute those points on the curve which lie to the positive side of the
% intersection point

display('slicing near coordinate axis intersection');

raw_slice_points = coordinate_slice(b_input, initial_start_points, ...
			random_line_coefficients, slice_value, variable_index, 'complex');

		
%find those paths which have already been tracked.
[points_connecting_to_a_current_intersection_point, raw_connections, reject_connecting_points] = reject_paths_which_go_to_unentered_intersection_point(...
	raw_slice_points, ...
	slice_value, intersection_points, ...
	b_input, variable_index);


connections.raw = raw_connections;


display('doing monodromy loops to determine unique paths');

%  go on a cauchy walk (monodromy loop) to determine which upper and lower
%  points connect to the intersection points.
[walked_paths, unconsumed_slice_points, unconsumed_connections] = cauchy_walk_complex(points_connecting_to_a_current_intersection_point, raw_connections, intersection_points, slice_value, component_degree, b_input, variable_index, num_sample_points, unique_point_threshold);
%upper_points and lower_points are returned here, having eliminated any
%points which were visited during the course of monodromy looping.


if ~isempty(unconsumed_slice_points)
	warning('there were unmapped points');
end

if ~isempty(unconsumed_connections)
	warning('there were unmapped connections');
end



slice_points.raw = raw_slice_points;

saveme.previous_paths = walked_paths;
saveme.previous_b_input = b_input_copyme;
saveme.previous_slice_value = slice_value; 
saveme.previous_slice_points = slice_points;
saveme.previous_connections = connections;
saveme.previous_reject_connecting_points = reject_connecting_points;  %#ok, saveme is indeed used.
save('previous_paths.mat','-struct','saveme')



%return to the parent directory.
cd ..
%leave the temp directory in place.


end







function [slice_points, connections, reject_connecting_points] = reject_paths_which_go_to_unentered_intersection_point(slice_points, slice_value, local_isect, b_input, variable_index)


connections = zeros(1,size(slice_points,2));

connecting_points = zeros(size(slice_points));

reject_connecting_points = [];

display('connecting slice points to intersection points');

for ii = 1:size(slice_points,2)
	a = connect_slice_to_intersection(b_input, slice_points(:,ii), slice_value, 0, variable_index,'complex');
	connecting_points(:,ii) = a;
end

%find the connecting point in the set of points we are interested in --
%local_isect.  
%if the connecting point is not found, then we should not do the rest of
%the calculations on the lower point.

delete_these_rows = [];
for ii = 1:size(connecting_points,2)
	ind = find_same_point(local_isect, connecting_points(:,ii), 1e-6, 1);
	% ind is nonempty if the point tracked to a point not in local_isect.
	if isempty(ind)
		delete_these_rows(end+1) = ii;
		reject_connecting_points(:,end+1) = connecting_points(:,ii);
		connections(ii) = 0;
	else
		connections(ii) = ind;
	end
end
slice_points(:,delete_these_rows) = [];
connections(delete_these_rows) = [];

reject_connecting_points = unique_up_to_threshold(reject_connecting_points,1e-8);

end
