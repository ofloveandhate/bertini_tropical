% parse a bertini1 input file into a struct.
% the struct potentially has fields:
%
%   config
%   variable_group
%   variable
%   hom_variable_group
%   function
%   subfunction
%   constant
%
%   where these fields have additional subfields.  the names of the
%   subfields are the names of objects appearing in the bertini input file.
%   the values of the subfields are their values, which may be empty.
%
%   for example, a system may have more than one variable_group.  in this
%   case, b_input.variable_group has two entries.  each of these entries
%   has subfields, the names of which are the names of the variables in
%   the variable_group.s  The values are the empty string ''.
%   variable_group x, y;
%   variable_group u, v;
%   
%    will have
%     
%             b_input.variable_group(1).x = '';
%             b_input.variable_group(1).y = '';
%             b_input.variable_group(2).u = '';
%             b_input.variable_group(2).v = '';
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

function [b_input] = parse(b_input, filename)

if nargin<2
	filename = 'input';
end

b_input.filename = filename;

file_status = exist(filename,'file');


if ~file_status
	error('filename ''%s'' is not an existing file',filename);
end

num_lines = getlinecount(filename);
fid = fopen(filename,'r');

line_number = 1; %initialize the counter for lines.

%the following indicators tell us whether we are inside the input or config
%sections, and help us tell whether an input file is valid.
finished_input = false;
finished_config = false;
started = false;
in_config = false;
in_input = true;

b_input.filename = filename;

if isempty(which('md5'))
	warning('missing the md5 function.  please add to the path a function ''md5'' which computes the md5 hash of a file, and returns that hash as a string.  This Bertini parser uses it to compare different input files without parsing.  Recommended function at https://www.mathworks.com/matlabcentral/fileexchange/5498-md5-signature-of-a-file');
else
	b_input.md5hash = md5(filename);
end
b_input.declared_symbols = {'I','special_number';'i','special_number';'Pi','special_number';'pi','special_number'};  

b_input.config = struct;
b_input.subfunction = cell(0,0);
last_line_of_file = NaN;

while ~(finished_input || line_number > num_lines)
	
	
	totalline = '';
	
	acceptable_input = false;
	
	while ~(acceptable_input || line_number > num_lines)
		currline = fgetl(fid);
		line_number = line_number+1;
		
		[locations] = strfind(currline,'%'); %find the beginning of a comment, if it exists
		
		if isempty(locations)
			ell = length(currline);
		else
			ell = locations(1)-1;
		end
		
		totalline = sprintf('%s%s',totalline,currline(1:ell)); %append the non-comment part of the string
		
		
		%trim off the leading spaces
		space_locations = isspace(totalline);
		if ~isempty(space_locations)
			if space_locations(1)
				totalline = totalline( find(space_locations~=1:length(space_locations), 1,'first' ) : end);
			end
		end


		%trim off the trailing spaces
		space_locations = isspace(totalline);  %creates a boolean array indicating whether characters are spaces
		while ~isempty(space_locations)
			if space_locations(end);
				totalline = totalline( 1: find(space_locations~=1:length(space_locations), 1,'last' )-1 );
			else
				break;
			end
			space_locations = isspace(totalline); %update.
		end
		
		
		if ~isempty(totalline)
			
			if strcmp(totalline,'CONFIG') || strcmp(totalline,'INPUT') ||strcmp(totalline(end),';')
				acceptable_input = true;
			else
% 				sprintf('''%s''',totalline(end))
			end

		end
		
		
	end % re: while
	
	
	if strcmp(totalline,'CONFIG')
		
		if in_input && started
			warning('beginning config while in INPUT section, around line %i',line_number);
			acceptable_input = false;
		end
		
		if finished_input || finished_config
			warning('beginning config after finished input or config, around line %i',line_number);
			acceptable_input = false;
		end
		
		if in_config
			warning('multiple CONFIG declarations, around line %i',line_number);
			acceptable_input = false;
		end
		in_config = true;
		in_input = false;
		
		
		
		continue
	end
	
	if strcmp(totalline,'INPUT')
		
		if in_config && ~finished_config
			warning('beginning INPUT before finished config section, around line %i',line_number);
			acceptable_input = false;
		end
		
		if ~in_config && ~finished_config
			finished_config = true;
		end
		
		if in_input
			warning('multiple INPUT declarations, around line %i',line_number);
			acceptable_input = false;
		end
		
		in_input = true;
		in_config = false;
		started = true;
		continue
	end
	
	
	
	
	
	if strcmp(totalline,'END;')
		if ~in_config && ~in_input
			error('stray END; around line %i',line_number);
		end
		
		if in_config
			finished_config = true;
			in_config = false;
			continue;
		end
		
		if in_input
			in_input = false;
			finished_input = true;
			last_line_of_file = line_number;
			continue;
		end
	end
		
	
	if in_config
		colon_location = strfind(totalline,':');
		
		if isempty(colon_location)
			error('configuration options must be separated from values by a colon');
		end
		
		option_name = totalline(1:colon_location-1);
		option_name = lower(option_name(~isspace(option_name)));
		
		option_value = totalline(colon_location+1:end-1);  %strip the colon and semicolon off.
		option_value = option_value(~isspace(option_value));
		
		b_input.config.(option_name) = option_value;
		continue;
	end
	
	if ~in_input
		continue;
	end
	
	if isempty(totalline)
		continue;
	end

	
	
	
	equals_locations = strfind(totalline,'=');
	
	
	
	if ~isempty(equals_locations)
		% there is an equals sign somewhere
		
		left_of_equals = totalline(1:equals_locations(1)-1);
		left_of_equals = left_of_equals(~isspace(left_of_equals));
		right_of_equals = totalline(equals_locations(1)+1:end);
		
		
		
		current_type = 'subfunction';
		for mm = 1: size(b_input.declared_symbols,1)
			if strcmp(left_of_equals, b_input.declared_symbols{mm,1})
				current_type = b_input.declared_symbols{mm,2};
				break;
			end
		end
		
		
		[start_indices,end_indices] = regexp(right_of_equals,'[a-zA-Z_][a-zA-Z0-9\[\]_]*\d*\.?\d+e[+-]?\d+');%[^a-zA-Z*][a-zA-Z_0-9*]
	

		undeclared_symbols_RHS = {};
		for jj = 1:length(start_indices)
			testme = right_of_equals(start_indices(jj):end_indices(jj));
			
			found_local = false;
			for mm = 1: size(b_input.declared_symbols,1)
				if strcmp(testme, b_input.declared_symbols{mm,1})
					found_local = true;
				end
			end
			
			if ~found_local
				undeclared_symbols_RHS{end+1} = testme;
			end
		end


		
		if ~isempty(undeclared_symbols_RHS)
			display(b_input.declared_symbols)
			display(undeclared_symbols_RHS);
			error('symbol inside ''%s'', on line %i',right_of_equals,line_number);
		end
		
		if strcmp(current_type, 'subfunction')
			[b_input] = declare_symbols(b_input, {left_of_equals}, 'subfunction');
		end
		
		index = 0;
		for qq = 1:size(b_input.(current_type),1)
			if strcmp(b_input.(current_type){qq,1},left_of_equals)
				index = qq;
				break;
			end
		end
		b_input.(current_type){index,1} = left_of_equals;
		b_input.(current_type){index,2} = totalline(equals_locations(1)+1:end-1);
		

		
		
	else  %empty equals_locations
		words = strsplit(totalline);
		first_word = words{1};
		
		remainder = totalline(length(first_word)+1:end-1);
		remainder = remainder(~isspace(remainder));
		
		words = strsplit(remainder,',');
		
		switch first_word
			case 'function'
				
			case 'variable_group'
				
			case 'hom_variable_group'
				
			case 'constant'
				
				
			otherwise
				error('declaration type %s not in lexicon',first_word);
		end
		
		[b_input] = declare_symbols(b_input, words, first_word);
		
		
		
		
		
	end
	


end


if isnan(last_line_of_file)
% 	warning('finished reading file, but INPUT never reached END;');
end




fclose(fid);



validate(b_input);





end












function [ok_flag,bad_symbols] = validate(b_input)

ok_flag = true;

bad_symbols = {};

for ii = 1:length(b_input.declared_symbols)
	n = b_input.declared_symbols{ii,1};
	t = b_input.declared_symbols{ii,2};

	if strcmp('special_number',t) || strcmp('variable_group',t) || strcmp('hom_variable_group',t) || strcmp('path_variable',t) || strcmp('variables',t)
		continue;
	end
	
	
	index = 0;
	for jj = 1:size(b_input.(t),1)
		if strcmp(b_input.(t){jj,1},n)
			index = jj;
			break;
		end
	end
	
	if isempty(b_input.(t){index,2})
		ok_flag = false;
		bad_symbols{end+1} = n;
	end
end


if ~ok_flag
	display('there were undefined symbols:');
	for ii = 1:length(bad_symbols)
		display(bad_symbols{ii});
	end
	
	error('cannot proceed');
end

end
