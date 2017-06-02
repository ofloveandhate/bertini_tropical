%  [] = write_bertini_input_file(variables, polynomial_system,varargin)
%
%  writes a bertini input file to the current folder.
%
%  arguments:

%    *  variables: you have three options.  1)a struct, whos fieldnames are the names of the
%    variables, and whos values are arbitrary.  2) a cell array of names of
%    variables. 3) a cell array of syms; that is, the variables themselves.

%    *  polynomial_system: a struct, whose names are the names of the
%    functions you want to write, and whose values are the strings of the
%    functions.

%    *  varargin:  an option-value pair, where the first element of the
%    pair is the name of the option, and whose value is the value of the
%    option.
%
%
%  options:
%    *  'filename': accompanied by a string, the value of which is the name
%    of the file you are creating.  'input' by default.

%    *  'options':  a struct, whose fieldnames are the names of the
%    bertini options you want to use, and whose values are the values of
%    those options. 

%    *  'constants': a struct, whose fieldnames are the names of constants,
%    and whose values are the values of the constants.

%    *  'parameters':  a struct containing the fieldname-value pairs of
%    parameters to be written.  use of parameters must be accompanied by
%    the use of the option 'userhomotopy' with non-zero value.

%    *  'subfunctions':  a struct containing the fieldname-value pairs of
%    subfunctions to be written.
%
%
%  the input file is written in the following order:
%    1) options into config section
%    2) pathvariable and parameter t and s declared, if userhomotopy.
%    3) variables declared, as variable_group if userhomotopy==0||2, or
%         variable if userhomotopy==1
%    4) constant
%    5) parameters
%    6) subfunctions
%    7) functions
%
%
%  noted shortfalls: cannot write homogeneous variable groups, nor multiple
%  variable groups.  only one variable group for right now.
%
% 
% syms x y
% variables = {x y};
% system.f = x^2 + y^2 - 1;
% options.tracktype = 1;
% write_bertini_input_file(variables,system,'options',options,'filename','exampleinput')

%  please inform the author immediately if this code doesn't do what you
%  think it should do.  
%
% copyright 2015, 2016, 2017 Dani Brake
% University of Notre Dame
% Applied and Computational Mathematics and Statistics
% University of Wisconsin
% Mathematics
% danielthebrake@gmail.com  brakeda@uwec.edu
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


function write_bertini_input_file(variables, polynomial_system,varargin)

options = struct();
filename='input';
constants = struct();
parameters = struct(); declare_parameters = false;
subfunctions = struct();
dependent_subfuncs = struct();

userhomotopy = 0;

for ii = 1:2:length(varargin)
	switch varargin{ii}
		case 'filename'
			filename = varargin{ii+1};
			
		case 'options'
			options = varargin{ii+1};
			
		case 'constants'
			constants = varargin{ii+1};
			
		case 'parameters'
			parameters = varargin{ii+1}; declare_parameters = true;
			
		case 'subfunctions'
			subfunctions = varargin{ii+1};
			
		case 'dependent_subfunctions'
			dependent_subfuncs = varargin{ii+1};
		otherwise 
			error('bad option ''%s'' to make_bertini_input_file_advanced',varargin{ii});
			
	end
end


%scan for use of the userhomotopy option.
f = fieldnames(options);
for ii = 1:length(f)
	if strcmpi(f{ii},'userhomotopy')
		if ischar(options.(f{ii}))
				userhomotopy = str2double(options.(f{ii}));
		else 
				userhomotopy = options.(f{ii});
		end
		
	end
end


%open the input file
fid = fopen(filename,'w');


% print the config section of the input file.
fprintf(fid,'CONFIG\n\n');
	print_definitions(fid,options,':');
fprintf(fid,'\n\nEND;\n\nINPUT\n\n');

%print the input section of the file.
if userhomotopy==0
	if isstruct(variables)
		print_definitions(fid,variables,'','variable_group');
	else
		print_declaration(fid,variables,'variable_group');
	end
elseif userhomotopy==1
	fprintf(fid,'pathvariable t;\nparameter s;\ns=t;\n\n');
	
	if isstruct(variables)
		print_definitions(fid,variables,'','variable');
	else
		print_declaration(fid,variables,'variable');
	end
	
else
	fprintf(fid,'pathvariable t;\nparameter s;\ns=t;\n\n');
	if isstruct(variables)
		print_definitions(fid,variables,'','variable_group');
	else
		print_declaration(fid,variables,'variable_group');
	end
end

%print the constants, parameteres, subfunctions, and system itself, in that
%order.
print_definitions(fid,constants,'=','constant');
if declare_parameters
	print_declaration(fid,parameters,'parameter');
end
print_definitions(fid,subfunctions,'=');

print_definitions(fid,dependent_subfuncs,'=');

print_definitions(fid,polynomial_system,'=','function');

%add the closing END statement
fprintf(fid,'\n\n\nEND;\n\n');

%add a time stamp
fprintf(fid,'%%this input file automatically generated\nby make_bertini_input_file\n%s\n\n\n',datestr(clock));

%finally close the input file
fclose(fid);


end



function print_definitions(fid,s,delimiter,name)

if nargin ==3
	if isstruct(s)
	print_struct_definitions(fid,s,delimiter);
	elseif iscell(s)
		print_cell_definitions(fid,s,delimiter);
	end
else
	if isstruct(s)
		print_struct_definitions(fid,s,delimiter,name);
	elseif iscell(s)
		print_cell_definitions(fid,s,delimiter,name);
	end
end

end




%  print a cell array to the file, declared as a 'name' type of object,
%  separated by delimiter.  s is the struct to write
function print_struct_definitions(fid,s,delimiter,name)

if ~isstruct(s)
	display(s);
	error('the above variable is not a struct...\ninput at print_struct for type %s is not a struct',name);
end

f = fieldnames(s);
if ~isempty(f)
	if nargin >=4
		
		fprintf(fid,'%s ',name);
		for ii = 1:length(f)
			fprintf(fid,'%s',f{ii});
			if ii == length(f)
				fprintf(fid,';\n');
			else
				fprintf(fid,', ');
			end
		end
		
	end
	
	if ~isempty(delimiter)
		for ii = 1:length(f)
			fprintf(fid,'%s %s ',f{ii},delimiter);
			print_value(fid,s.(f{ii}));
			fprintf(fid,';\n');
		end
	end
	
	fprintf(fid,'\n');
end


end


function print_cell_definitions(fid,s,delimiter,name)

if ~iscell(s)
	display(s);
	error('the above variable is not a cell...\ninput at print_cell for type %s is not a cell',name);
end

num_def = size(s,1);
if ~isempty(s)
	if nargin >=4
		
		
		fprintf(fid,'%s ',name);
		for ii = 1:num_def
			fprintf(fid,'%s',s{ii,1});
			if ii == num_def
				fprintf(fid,';\n');
			else
				fprintf(fid,', ');
			end
		end
		
	end
	
	if ~isempty(delimiter)
		for ii = 1:num_def
			fprintf(fid,'%s %s ',s{ii,1},delimiter);
			print_value(fid,s{ii,2});
			fprintf(fid,';\n');
		end
	end
	
	fprintf(fid,'\n');
end


end





%  print a cell array to the file, declared as a 'name' type of object.
function print_declaration(fid,arr,name)
fprintf(fid,'%s ',name);
for ii = 1:length(arr)
	print_value(fid,arr{ii});
	if ii == length(arr)
		fprintf(fid,';\n');
	else
		fprintf(fid,', ');
	end
end

end


function print_value(fid,value)

if ischar(value)
	fprintf(fid,'%s',value);
elseif isint(value)
	fprintf(fid,'%i',value);
elseif isfloat(value)
	if isreal(value)
		fprintf(fid,'%1.16f',value);
	else
		fprintf(fid,'%1.16f+I*%1.16f',real(value),imag(value));
	end
elseif isa(value,'sym')
	if isreal(value)
		fprintf(fid,'%s',char(value));
	else

		foundvars = symvar(value);
		
		if length(foundvars)>=1
			
			if and(length(foundvars)==1,strcmp(char(foundvars(1)),'i'))
				fprintf(fid,'%s+I*%s',char(real(value)),char(imag(value)));
			else
				fprintf(fid,'%s',char(value));
			end
		else
			fprintf(fid,'%s+I*%s',char(real(value)),char(imag(value)));
		end
		
	end
else
	error('unknown type of object')
end
end
