%parse from a witness_data file from Bertini into a witness_data class
%object defined in Matlab
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

function data = parse(data, filename)

if nargin==1
	filename = 'witness_data';
end


fid = fopen(filename,'r');


reset(data);

num_vars = fscanf(fid,'%i',[1 1]);
data.num_vars = num_vars;

num_nonempty_codims = fscanf(fid,'%i',[1 1]);

data.nonempty_codimensions = zeros(1,num_nonempty_codims);

for ii = 1:num_nonempty_codims
	
	codim = fscanf(fid,'%i',[1 1]);
	data.nonempty_codimensions(ii) = codim;
	
	num_pts = fscanf(fid,'%i',[1 1]);
	
	
	components_encountered = [];

	current_witness_sets = repmat(witness_set,[0 0]);
	
	for zz = 1:num_pts
		precision = fscanf(fid,'%i\n',[1 1]);
		point = vpa(zeros(num_vars,1));
		for jj = 1:num_vars
			tmp1 = vpa(fscanf(fid,'%s ',[1 1]));
			tmp2 = vpa(fscanf(fid,'%s\n',[1 1]));
 			point(jj) = tmp1+1i*tmp2;
		end
		precision_last_approx = fscanf(fid,'%i',[1 1]);
		prev_point = zeros(num_vars,1);
		for jj = 1:num_vars
			tmp = fscanf(fid,'%f',[2 1]);
			prev_point(jj) = tmp(1)+1i*tmp(2);
		end

		condition_number = fscanf(fid,'%f',[1 1]);
		corank_of_jacobian = fscanf(fid,'%i',[1 1]);
		smallest_nonzero_sing_val = fscanf(fid,'%f',[1 1]);
		largest_zero_sing_val = fscanf(fid,'%f',[1 1]);

		type = fscanf(fid,'%i',[1 1]); %10==nonsingular, 15==singular

		multiplicity = fscanf(fid,'%i',[1 1]);
		component_number = fscanf(fid,'%i',[1 1]);
		deflations_needed = fscanf(fid,'%i',[1 1]);

		if component_number<0
			display('detected unclassified component')
				continue
		end
			
		if isempty(find(components_encountered==component_number, 1))
			
			components_encountered(end+1) = component_number; %store the fact that we encountered this component number
			
			current_witness_sets(component_number+1) = witness_set; % create empty witness set
			current_witness_sets(component_number+1).component_number = component_number; %assign it this number
			current_witness_sets(component_number+1).codimension = codim;
			
			if type==10
				current_witness_sets(component_number+1).singular = false;
			elseif type==15
				current_witness_sets(component_number+1).singular = true;
			else
				error('unexpected ''type''');
			end
			
		else
			
		end
		
		current_witness_sets(component_number+1).add_point(point);

	end
	
	if isempty(data.witness_sets)
		data.witness_sets = current_witness_sets;
	else
		data.witness_sets = [data.witness_sets current_witness_sets];
	end
end

parity_check = fscanf(fid,'%i', [1 1]);


if parity_check~=-1
	error('did not get -1 after the point block');
end



num_type = fscanf(fid,'%i',[1 1]);

for ii = 1:num_nonempty_codims

	codim = data.nonempty_codimensions(ii);
	% the codimensions are expected to  come in the same order as above.

	%first, get the randomization matrix
	
	
	num_rows = fscanf(fid,'%i',[1 1]);
	num_cols = fscanf(fid,'%i',[1 1]);
	
	if num_type==2
		randomization_matrix = vpa(zeros(num_rows, num_cols));
	else
		randomization_matrix = zeros(num_rows, num_cols);
	end
	for jj = 1:num_rows
		for kk = 1:num_cols
			if num_type==2
				tmp1 = vpa(fscanf(fid,'%s',[1 1]));
				tmp2 = vpa(fscanf(fid,'%s',[1 1]));
			else
				tmp1 = fscanf(fid,'%f',[1 1]);
				tmp2 = fscanf(fid,'%f',[1 1]);
			end
			
			randomization_matrix(jj,kk) = tmp1+1i*tmp2;
			
		end
	end
	
	for zz = 1:length(data.witness_sets)
		if data.witness_sets(zz).codimension == codim
			data.witness_sets(zz).set_randomization_matrix(randomization_matrix);
		end
	end
	
	if num_type==2
		hom_matrix = vpa(zeros(num_rows, num_cols));
	else
		hom_matrix = zeros(num_rows, num_cols);
	end
	
	%get hom matrix
	for jj = 1:num_rows
		for kk = 1:num_cols
			% according to the book, these entries are all integers...
			if num_type==2
				tmp1 = vpa(fscanf(fid,'%s',[1 1]));
			else
				tmp1 = fscanf(fid,'%s',[1 1]);
			end
			hom_matrix(jj,kk) = tmp1;
		end
	end
	
	
	% vector h for homogenization
	num_entries_in_hom_vector = fscanf(fid,'%i',[1 1]);
	
	
	if num_type==2
		hom_vector = vpa(zeros(num_entries_in_hom_vector,1));
	else
		hom_vector = zeros(num_entries_in_hom_vector,1);
	end
	
	for jj = 1:num_entries_in_hom_vector
		if num_type==2
			tmp1 = vpa(fscanf(fid,'%s',[1 1]));
			tmp2 = vpa(fscanf(fid,'%s',[1 1]));
		else
			tmp1 = vpa(fscanf(fid,'%s',[1 1]));
			tmp2 = vpa(fscanf(fid,'%s',[1 1]));
		end
		hom_vector(jj) = tmp1+1i*tmp2;
	end
	
	
	% get the hom_var_const...
	tmp1 = vpa(fscanf(fid,'%s',[1 1]));
	tmp2 = vpa(fscanf(fid,'%s',[1 1]));
	
	hom_var_const = tmp1 + 1i*tmp2;
	%according to the book, is 0 if affine, random complex if projective
	
	
	%get the linear slices
	num_rows = fscanf(fid,'%i',[1 1]); % the number of slices there are
	num_cols = fscanf(fid,'%i',[1 1]); % the number of variables there are
	
	for jj = 1:num_rows
		if num_type==2
			linear = vpa(zeros(num_vars,1));
		else
			linear = zeros(num_vars,1);
		end
		
		for kk = 1:num_cols
			if num_type==2
				tmp1 = vpa(fscanf(fid,'%s',[1 1]));
				tmp2 = vpa(fscanf(fid,'%s',[1 1]));
			else
				tmp1 = fscanf(fid,'%f',[1 1]);
				tmp2 = fscanf(fid,'%f',[1 1]);
			end
			linear(kk) = tmp1+1i*tmp2;
		end
		
		for zz = 1:length(data.witness_sets)
			if data.witness_sets(zz).codimension == codim
				data.witness_sets(zz).add_linear(linear);
			end
		end
	end
	
	
	size_patch = fscanf(fid,'%i', [1 1]);
	
	if num_type==2
		patch = vpa(zeros(size_patch,1));
	else
		patch = zeros(size_patch,1);
	end
	for jj = 1:size_patch
		if num_type==2
			tmp1 = fscanf(fid,'%s',[1 1]);
			tmp2 = fscanf(fid,'%s',[1 1]);

			s1 = find(tmp1=='/',1);
			n1 = vpa(tmp1(1:s1-1));
			d1 = vpa(tmp1(s1+1:end));

			s2 = find(tmp2=='/',1);
			n2 = vpa(tmp2(1:s2-1));
			d2 = vpa(tmp2(s2+1:end));

			patch(jj) = n1/d1+1i*n2/d2;
		else
			tmp3 = fscanf(fid,'%f',[1 1]);
			tmp4 = fscanf(fid,'%f',[1 1]);
			patch(jj) = tmp3+1i*tmp4;
		end
	end
	
	for zz = 1:length(data.witness_sets)
		if data.witness_sets(zz).codimension == codim
			data.witness_sets(zz).add_patch(patch);
		end
	end
	
end




fclose(fid);


end
