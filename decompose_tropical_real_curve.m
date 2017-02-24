%%  [produced_rays, ray_multiplicities, aux_data] = tropical_real(varargin)
%
%  a function which decomposes a real curve's tropicalization.
%
%  the input for this function is a bertini input file, whos name is input
%  by default.
%
%  the main two output variables are produced_rays, and
%  ray_multiplicities.  all other data is output via the aux_data
%  struct, and includes cauchy paths, puiseux series coefficient,
%  and a number of other things
%
%  options are passed to this function by string-valued option names, as
%  detailed below.
%
% command-line options:
%
% * 'input' - a string specifying the name of the Bertini input file to use
% * 'options' - a struct, the fieldnames of which are bertini options, and the
% * values of which are the values of those options
% * 'purge' - boolean.  true means delete previously existing temporary
% directories prior to starting, false is leave them in place.  false is
% default.
% * 'defaultslicevalue' - positive float, the default value for slicing near
% coordinate axes.  the slice value used is the lesser of (this value, and
% the nearest absolute value of critical points).  the default for this
% setting was arbitrarily chosen to be 0.1.
% * 'numsamplepoints' - even integer, determines the number of samples
% taken for monodromy around the coordinate axes.  default is 8.  larger
% problems may require a higher number.
% * 'puiseuxthreshold' - computed puiseux series coefficients less than this
% tolerance are considered zero.
% * 'intersectionthreshold' - coordinates of intersection points less than this
% tolerance are considered zero.
% * 'intersectionuniquethreshold' - intersection points are unique if
% separated by this distance or more.
% * 'cauchyuniquethreshold' - points on cauchy loops are unique if separated
% by more than this distance.
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



function [produced_rays, ray_multiplicities, aux_data] = decompose_tropical_real_curve(varargin)





%%%%%%%%%%%%%%%%%%%%%
%
% preliminary setup
%
%%%%%%%%%%%%%%%%%%%%%


%% parse the command line and read the input file
[curve_system_name, cauchy_walk_point_uniqueness_threshold, intersection_point_uniqueness_threshold, intersection_zero_coord_threshold, puiseux_zero_coeff_threshold,min_slice_value,default_slice_value,remove_prev_dir,num_sample_points,user_config] = parse_varargin(varargin{:});



b_input = bertini_input(curve_system_name);
if isfield(b_input.config,'tracktype')
	b_input.config = rmfield(b_input.config,'tracktype');
end
f = fieldnames(user_config);
for ii = 1:length(f)
	b_input.config.( lower(f{ii}) ) = user_config.(f{ii});
end

system_variables = b_input.variable_group;
num_vars = numel(b_input.variable_group);





%% get the numerical irreducible decomposition of the system
if isempty(dir('witness_data'))
	error('cannot continue without bertini file ''witness_data''.  please generate it.');
end

w_data = witness_data('parse','dehom');

selections = w_data.unidimensionsal_select(1); % choose only from one-dimensional components
w = w_data.construct_agglomeration(selections);

aux_data.initial_start_points = w.points;
aux_data.random_line_coefficients = w.linears;

component_degree = size(aux_data.initial_start_points,2); %the degree of the complex component

starting_directory = move_into_working_directory(remove_prev_dir);



%%%%%%%%%%%%%%
%
%  now starts the actual work
%
%%%%%%%%%%%%%%%%%%





%% get only the unique intersection points.

[raw_intersection_points] = compute_intersection_points(b_input, aux_data,'real');


%this is here, and not in the compute_intersection_points function, so that
%the intersections points can be thresholded differently if the setting
%changes on a subsequent call.
[intersection_points,~] = unique_up_to_threshold(raw_intersection_points, intersection_point_uniqueness_threshold); % ~ replaced indices


%set up a bunch of empty variables which grow during the
%variable-by-variable loop which follows.
raw_corresponding_coefficients = {};
raw_corresponding_variables = {};
raw_rays = [];
considered_points = [];
aux_data.paths = [];





% compute the critical points with respect to the current time variable
[aux_data.generic_critical_start_points, detjac_input] = compute_critical_points_start_points(b_input,aux_data.initial_start_points,aux_data.random_line_coefficients);



%cache a copy into isect for use and modification throughout the next loop.
isect = intersection_points;

slice_values = zeros(num_vars,1);


%% we do a variable-by-variable study, computing the rays for each intersection point.

for ii= 1:num_vars
	
	%% get the indices of those intersection points such that variable ii is 0
	intersect_this_variable = find(abs(isect(ii,:))<intersection_zero_coord_threshold);
	
	
	% put the found points into a matrix for use in this iteration of the
	% loop.  then delete from the isect variable, so that we don't study the
	% intersection points more than once.
	local_isect = isect(:,intersect_this_variable);
	isect(:,intersect_this_variable) = []; %delete the used points.  don't worry, we have a cached copy up top.
	
	
	% if there are no intersection points for this variable, then there's
	% nothing to do.  continue.
	if isempty(intersect_this_variable)
		fprintf('\n\n\nno unstudied intersection points for variable %s\n\n\n\n',system_variables{ii});
		slice_values(ii) = nan;
		continue;
	end
	
	
	
	%% compute the critical points with respect to the current time variable
	aux_data.finite_critical_points.(system_variables{ii}) = compute_specific_critical_points(detjac_input,ii,aux_data.generic_critical_start_points);
	
	
	% by default slice at default (arbitrary, but up to user to override) place.
	slice_value = default_slice_value;
	
	if ~isempty(aux_data.finite_critical_points.(system_variables{ii}))
		%the following syntax is likely broken on older versions of matlab.
		crit_point_values = aux_data.finite_critical_points.(system_variables{ii})(ii,:); % project onto the kth variable
		
		%threshold away critical values at 0.
		crit_point_values = crit_point_values( abs(crit_point_values)>intersection_zero_coord_threshold );
		
		nonzero_intersection_coordinates = real(intersection_points(ii,:));
		nonzero_intersection_coordinates = nonzero_intersection_coordinates(abs(nonzero_intersection_coordinates)>intersection_zero_coord_threshold);
		
		if ~isempty(crit_point_values)
			slice_value = max(min([default_slice_value abs(crit_point_values)/2 abs(nonzero_intersection_coordinates)/2]),min_slice_value);
		end
		
	end
	
	slice_values(ii) = slice_value;
	fprintf('slicing value for t=%s is %1.4e\n',system_variables{ii},slice_value)
	
	
	%% slice the curve near the intersection points, and then do monodromy
	%loops around the intersection points to determine the number of unique
	%paths leading to the point.
	%
	%this call builds the paths, and stores them.  to the paths are
	%attached a bunch of data -- the cycle number, the center point, the previously known
	%points which were visited, and the path itself.
	[tmp_paths,rejected_path_intersection_points,slice_points, connections] = slice_and_monodromy_real(slice_value, local_isect, aux_data.initial_start_points, aux_data.random_line_coefficients, component_degree, b_input, ii, num_sample_points, cauchy_walk_point_uniqueness_threshold);
	
	if size(unique_up_to_threshold([considered_points rejected_path_intersection_points],intersection_point_uniqueness_threshold),2) ~= size(considered_points,2)
		warning('have encountered intersection points which were not previously studied.  the results of this decomposition are almost certainly incorrect.  these points have been stored.');
		aux_data.rejected_path_intersection_points.(system_variables{ii}) = rejected_path_intersection_points;
	end
	
	
	%add the paths which were computed.
	aux_data.paths = [aux_data.paths tmp_paths];
	aux_data.connections.(system_variables{ii}) = connections;
	aux_data.studied_intersection_points.(system_variables{ii}) = local_isect;
	aux_data.slice_points.(system_variables{ii}) = slice_points;
	
	
	
	%preallocate
	numerators = zeros(num_vars,1);
	
	%%  compute the rays using cauchy integrals
	% for each path walked during slice and separate.
	for zz = 1:length(tmp_paths)
		computed_coefficients = zeros(num_vars,1);
		
		%grab the current center point.
		curr_center_point = tmp_paths(zz).center_point;
		
		%now to compute the cauchy integrals, for the puiseux coefficients.
		for vv=1:num_vars
			
			if 0 %abs(curr_center_point(vv))>intersection_zero_coord_threshold
				%a non-zero coordinate in the center point of a loop
				%indicates that the puiseux series has a constant term.
				%this is k=0, and in this case, the puiseux coefficient
				%is the center point coordinate.
				k = 0;
				computed_coeff = curr_center_point(vv);
			else
				
				
				
				k = 0; %initialize
				computed_coeff = cauchy_integral(tmp_paths(zz).path(vv,:), tmp_paths(zz).radius, tmp_paths(zz).cycle_number, k);
				
				largest_coeff = [computed_coeff k];
				
				while and(abs(computed_coeff)<puiseux_zero_coeff_threshold, k<=component_degree)
					k = k+1;
					computed_coeff = cauchy_integral(tmp_paths(zz).path(vv,:), tmp_paths(zz).radius, tmp_paths(zz).cycle_number, k);
					if abs(computed_coeff) > abs(largest_coeff(1))
						largest_coeff = [computed_coeff k];
					end
				end
				
				if k> component_degree
					warning('computed valuation for variable %i exceeds maximum theoretical valuation. this is certainly an error. largest puiseux coefficient was %e, for term %i press any key to continue.',vv,largest_coeff(1),largest_coeff(2));
					pause
				end
			end % re: if abs()
			
			
			computed_coefficients(vv) = computed_coeff;
			numerators(vv) = k;
		end
		
		raw_rays(1:num_vars,end+1) = -numerators;

		raw_corresponding_variables{end+1} = ii;
		raw_corresponding_coefficients{end+1} = computed_coefficients;
	end
	
	considered_points = [considered_points local_isect];
end

if size(considered_points,2) ~= size(intersection_points,2)
	display('the number of considered points does not match the total number of intersection points!');
	pause;
end


cd(starting_directory);


%% compute the number of times each unique ray appears in raw_rays.

produced_rays = unique_up_to_threshold(raw_rays, 1e-4);
ray_multiplicities = zeros(1,size(produced_rays,2));
corresponding_variables = cell(1,size(produced_rays,2));
corresponding_coefficients = cell(1,size(produced_rays,2));
corresponding_paths = cell(1,size(produced_rays,2));
corresponding_centers = cell(1,size(produced_rays,2));


for ii = 1:size(produced_rays,2)
	
	%first, copy the current unique ray, and repmat it.  then subtract
	%from raw_rays, and count the number of coordinates of that
	%difference which were 0.
	a = sum(abs(repmat(produced_rays(:,ii) , [1 size(raw_rays,2)]) - raw_rays) < 1e-4,1);
	
	
	% this gets those indices of raw_rays (columns of) which were
	% exactly the current (ii) produced_ray.
	% todo: replace this call with an all() call
	ind = find(a==num_vars);
	
	for jj = 1:length(ind)
		corresponding_variables{ii} = [corresponding_variables{ii} raw_corresponding_variables{ind(jj)}];
		corresponding_coefficients{ii} = [corresponding_coefficients{ii} raw_corresponding_coefficients(ind(jj))];
		corresponding_centers{ii} = [corresponding_centers{ii} aux_data.paths(ind(jj)).center_point];
		corresponding_paths{ii} = ind; % which paths led to the computation of the iith produced_ray.
	end
	
	g = gcd_all(produced_rays(:,ii)); %reduce
	produced_rays(:,ii) = produced_rays(:,ii)/g;
	ray_multiplicities(ii) = g * sum( a==num_vars );
end



if ~isempty(aux_data.paths)
	centers = [aux_data.paths(:).center_point];
	
	for ii = 1:size(intersection_points,2)
		if isempty(find_same_point(centers,intersection_points(:,ii), intersection_point_uniqueness_threshold,1))
			warning('intersection point %i was not the center of any path',ii);
		end
	end
end


aux_data.raw_intersection_points = raw_intersection_points;
aux_data.slice_values = slice_values;
aux_data.raw_rays = raw_rays;
aux_data.intersection_points = intersection_points;
aux_data.corresponding_variables = corresponding_variables;
aux_data.corresponding_coefficients = corresponding_coefficients;
aux_data.corresponding_paths = corresponding_paths;
aux_data.corresponding_centers = corresponding_centers;
end







%% change directories into a working directory

function starting_directory = move_into_working_directory(remove_prev_dir)

starting_directory = pwd;
[~, folder_name, ~] = fileparts(starting_directory);

% set up and and move into a temporary directory, so we keep the directory
% relatively clean.  otherwise there's an explosion of files.
dirname = sprintf('%s_real_tropical_decomposition',folder_name);

if remove_prev_dir
	dirlist = dir([dirname '*']);
	for ii = 1:length(dirlist)
		rmdir(dirlist(ii).name,'s')
	end
end

if ~isdir(dirname)
	mkdir(dirname);
end
cd(dirname);


end



%%  set program parameters to defaults, and parse the command line options

function [curve_system_name, cauchy_walk_point_uniqueness_threshold, intersection_point_uniqueness_threshold, intersection_zero_coord_threshold, puiseux_zero_coeff_threshold,min_slice_value,default_slice_value,remove_prev_dir,num_sample_points,user_config] = parse_varargin(varargin)


curve_system_name = 'input';

intersection_point_uniqueness_threshold = 1e-8;
cauchy_walk_point_uniqueness_threshold = 1e-8;
intersection_zero_coord_threshold = 1e-8;
puiseux_zero_coeff_threshold = 1e-10;

min_slice_value = 1e-5;
default_slice_value = 0.1;

remove_prev_dir = false;

num_sample_points = 16; %this number MUST be even

user_config = struct;

option_counter = 0;

%% parse out the options from varargin
while option_counter < length(varargin)
	option_counter = option_counter+1;
	curr_opt_name = varargin{option_counter};
	
	switch curr_opt_name
		case 'input'
			option_counter = option_counter+1;
			curve_system_name = varargin{option_counter};
			
		case 'options'
			option_counter = option_counter+1;
			user_config = varargin{option_counter};
			if isfield(user_config,'tracktype')
				error('do not specify a track type');
			end
		case 'purge'
			remove_prev_dir = true;
			
		case 'defaultslicevalue'
			option_counter = option_counter+1;
			curve_system_name = varargin{option_counter};
			
		case 'numsamplepoints'
			option_counter = option_counter+1;
			num_sample_points = varargin{option_counter};
			if mod(num_sample_points,2)~=0
				error('the number of sample points for monodromy MUST be even');
			end
			
		case 'puiseuxthreshold'
			option_counter = option_counter+1;
			
			puiseux_zero_coeff_threshold = varargin{option_counter};
			
			if ~isnumeric( puiseux_zero_coeff_threshold )
				error('threshold for puiseux coefficients must be a number');
			end
			
			if puiseux_zero_coeff_threshold < 0
				error('zero-threshold for puiseux coefficients must be greater than 0.')
			end
			
		case 'intersectionthreshold'
			option_counter = option_counter+1;
			
			intersection_zero_coord_threshold = varargin{option_counter};
			
			if ~isnumeric( intersection_zero_coord_threshold )
				error('threshold for intersection coordinates must be a number');
			end
			
			if intersection_zero_coord_threshold < 0
				error('zero-threshold for intersection coordinates must be greater than 0.')
			end
			
		case 'intersectionuniquethreshold'
			option_counter = option_counter+1;
			
			intersection_point_uniqueness_threshold = varargin{option_counter};
			
			if ~isnumeric( intersection_point_uniqueness_threshold )
				error('uniqueness threshold for intersection points must be a number');
			end
			
			if intersection_point_uniqueness_threshold < 0
				error('uniqueness threshold for intersection points must be greater than 0.')
			end
			
			
		case 'cauchyuniquethreshold'
			option_counter = option_counter+1;
			
			cauchy_walk_point_uniqueness_threshold = varargin{option_counter};
			
			if ~isnumeric( cauchy_walk_point_uniqueness_threshold )
				error('uniqueness threshold for points on cauchy loop must be a number');
			end
			
			if cauchy_walk_point_uniqueness_threshold < 0
				error('uniqueness threshold for points on cauchy loop must be greater than 0.')
			end
			
			
		case 'minslicevalue'
			option_counter = option_counter+1;
			
			min_slice_value = varargin{option_counter};
			
			if ~isnumeric( min_slice_value )
				error('minimum slice value must be a number');
			end
			
			if min_slice_value < 0
				error('minimum slice value must be greater than 0.')
			end
			
			
		otherwise
			error('bad option name %s',curr_opt_name);
	end
end


end



