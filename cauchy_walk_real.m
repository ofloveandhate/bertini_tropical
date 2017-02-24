%  walk around a pole in complex time.  the pole is where one variable
%  takes on value 0.  the purpose of the walk is to discover which other
%  points lie on the same curve around the pole.
%
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


function [walked_paths, unvisited_starting_points, unvisited_starting_connections, unvisited_antipodal_points, unvisited_antipodal_connections] = cauchy_walk_real(starting_points, starting_connections, antipodal_points, antipodal_connections, center_points, slice_value, component_degree, b_input_copyme, variable_index, num_sample_points, unique_point_threshold)

thresh = unique_point_threshold;


b_input = copy(b_input_copyme);


b_input.config.UserHomotopy = 2;

system_variables = b_input.variable_group;
num_vars = numel(system_variables);

if mod(num_sample_points,2)~=0
	error('num_sample_points MUST be even, because we have to walk half way around the pole');
end

b_input.declare_symbols({'walk_around'},'function');
b_input.declare_and_define('eps',slice_value,'constant');
	
	
num_paths = 0;

walked_paths = repmat(cauchy_path(),1,size(starting_points,2));

unvisited_starting_connections = starting_connections;
unvisited_antipodal_connections = antipodal_connections;

unvisited_starting_points = starting_points; % take a copy
unvisited_antipodal_points = antipodal_points;

while size(unvisited_starting_points,2)>0 %we delete points from unvisited_starting_points as this loop proceeds.

	num_paths = num_paths+1;
	walked_paths(num_paths) = cauchy_path();  % initialize to new class object.
	walked_path = walked_paths(num_paths);  %using handle class here, so that modification of one is propagated to the other.
	walked_path.path_variable = variable_index; %should this be the name or symbol rather than index?
	
	
	walked_path.set_start_point(unvisited_starting_points(:,1));
	walked_path.center_point = center_points(:,unvisited_starting_connections(1));
	walked_path.set_radius(abs(slice_value));
	%walk around t = eps to t = eps in complex time
	current_point = unvisited_starting_points(:,1);
	periodic_point = current_point;
	
	
	%initialize some data containers
	visited_starting_points = [];
	visited_antipodal_points = current_point;
	walked_path.add_connecting_point(current_point);
	
	unvisited_starting_connections(1) = [];
	unvisited_starting_points(:,1) = []; %delete the point from the list
	
	cycle_number = 0;
	iteration = 0;
	
	%enter the loop for walking around t=0 in complex time
	while iteration<component_degree
		
		iteration = iteration+1; % increment this counter
		
		
		for jj = 1:num_sample_points
			
			startname = sprintf('start%i_%i',iteration,jj);
			inputname = sprintf('input_walk_%i_%i',iteration,jj);
			
			
			b_input.define_symbol('walk_around',...
				sprintf('%s - eps*exp( I*(s*2*Pi*%i/%i +  (1-s)*2*Pi*%i/%i) )',system_variables{variable_index}, jj-1,num_sample_points,jj,num_sample_points));
			
			
			write_generic_solns(current_point,startname);
			
			write_bertini_input_file(system_variables, b_input.functions,'filename',inputname,'options',b_input.config,'constants',b_input.constant,'subfunctions',b_input.subfunction);%
			
			
			bertini('filename',inputname,'startname',startname,'stifle','rerunonfail',10); %,'mpi',min(4,size(current_point,2))
			
			next_point = get_generic_solns('nonsingular_solutions',num_vars);
			
			movefile('nonsingular_solutions',['nonsingular_solutions_' num2str(iteration) '_' num2str(jj)])
			
			if isempty(next_point)
				next_point = get_generic_solns('finite_solutions',num_vars);
				if ~isempty(next_point)
					warning('empty point in cauchy path...read from finite_solutions instead of nonsingular_solutions.consider changing condnumthreshold.');
				else
					error('path failure in cauchy_path tracking.  consult a Bertini expert.');
				end
				
			end
			
			walked_path.add_point_to_path(next_point);
			
			
			current_point = next_point;
			
			%halfway around the loop, we check and see if we visited any
			%antipodal points.
			if jj==num_sample_points/2
				visited_index = find_same_point(unvisited_antipodal_points, current_point, unique_point_threshold);
				 
				if isempty(visited_index)
					continue
				end
				
				if length(visited_index)>1
					thresh = unique_point_threshold;
					while length(visited_index)>1
						thresh = thresh/2;
						visited_index = find_same_point(unvisited_antipodal_points, current_point, thresh);

						if isempty(visited_index)
							thresh = thresh*sqrt(5);
						end
					end
					
					warning('had to shrink unique_point_threshold to %1.4e to get a unique visited antipodal point on monodromy loop',thresh);
					
				end
		
				walked_path.add_connecting_point(unvisited_antipodal_points(:,visited_index));
				visited_antipodal_points = [visited_antipodal_points unvisited_antipodal_points(:,visited_index)];
				
				unvisited_antipodal_connections = unvisited_antipodal_connections(1:end~=visited_index);
				unvisited_antipodal_points = unvisited_antipodal_points(:,1:end~=visited_index); %delete the visited point.
				
			end
		end %re: for jj
		
		cycle_number = cycle_number+1;
		
		%now we check to see if we came back to where we started.  if so,
		%then done with this particular loop, and will either be done, or
		%start a new loop on the first unvisited starting point if it
		%exists.
		if ~isempty(find_same_point(periodic_point, current_point, unique_point_threshold))
			break;
		end
		
		
		%got here, so not back at the beginning of the loop.
		
		
		visited_index = find_same_point(unvisited_starting_points, current_point, unique_point_threshold);
		
		
		%are we at a previously encountered point?  if not, we skip the
		%below things.
		if isempty(visited_index)
			continue
		end
		
		
		%well, we're here, so did not return to the start of the path, and
		%are at one of the points on the same line as the starting point.
		
		if length(visited_index)>1
			thresh = unique_point_threshold;
			while length(visited_index)>1
				thresh = thresh/2;
				visited_index = find_same_point(unvisited_starting_points, current_point, thresh);
				
				if isempty(visited_index)
					thresh = thresh*sqrt(5);
				end
			end
			
			warning('had to shrink unique_point_threshold to %1.4e\nto get a unique visited starting point on monodromy loop\n',thresh);
		end
		

		walked_path.add_connecting_point(unvisited_starting_points(:,visited_index));
		visited_starting_points = [visited_starting_points unvisited_starting_points(:,visited_index)];
		unvisited_starting_points = unvisited_starting_points(:,1:end~=visited_index);  %delete the visited points.
		unvisited_starting_connections = unvisited_starting_connections(1:end~=visited_index);
		
		
% 		unvisited_antipodal_points
% 		unvisited_starting_points
% 		visited_starting_points
% 		visited_antipodal_points
% 		max(abs(periodic_point-current_point))
% 		

		
	end %re: while and(max(...
	
	%need to do better error catching here for challenging decompositions.
	if isempty(find_same_point(periodic_point, current_point, thresh))
		warning('exceeded maximum theoretical upper bound on cycle number\ndid not return to start in %i cycles',component_degree)
		pause
	end
	
	walked_path.cycle_number = cycle_number;
	
end %re: while size(unvisited_starting_points...

walked_paths = walked_paths(1:num_paths);
end %re: function cauchy_walk














