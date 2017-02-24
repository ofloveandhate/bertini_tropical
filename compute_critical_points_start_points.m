function [generic_start_points, detjac_input] = compute_critical_points_start_points(b_input_copyme, initial_points, random_line_coefficients)
%
%
% compute generic critical points with respect to a particular variable, for a
% tropical curve.
%
%
% if successful, a .mat file named successful_run.mat is created, having
% the following data:
% 
% previous_b_input = b_input_copyme
% previous_generic_critical_points = generic_critical_points
% previous_input_for_moving = input_for_moving_to_particular_coefficients
% previous_detjac = detjac_input
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

dirname = sprintf('critpts_buildup');

if exist(dirname,'dir')==0
	mkdir(dirname)
end
cd(dirname)




if exist('successful_run.mat','file')
	load('successful_run.mat');
	
	if issame_system(previous_b_input,b_input_copyme)
		display('using previously generated generic critpt calculation');
		generic_start_points = previous_generic_start_points;
		detjac_input = previous_detjac_input;
		cd ..
		return
	end
end


b_input = copy(b_input_copyme);

system_variables = b_input.variable_group;
num_vars = length(system_variables);



linmove_input = b_input;  %cache a copy to be used for the linear moves




function_degrees = [get_function_degrees(copy(b_input));1];









%detjac_input gets some constants and subfunctions added to it.
%
%  constant.DETJAC_COEFF_varname
%  subfunction.(sprintf('DIFF_%s_%s'))

[detjac_input] = compute_detjac_function(copy(b_input));








display('regenerating to detjac function for generic critical points');

%create an initial empty witness set to contain the regeneration witness
%points.
W = witness_set;

M = [rand_complex(num_vars,num_vars); rand_complex(1,num_vars)];  %_complex random matrix for inversion for the nullvector start values.  the last row is the patch.
for ii = 1:num_vars+1
	for jj = 1:num_vars
		detjac_input.declare_and_define(sprintf('M_%i_%i',ii,jj),M(ii,jj),'constant');
	end
end


b = zeros(num_vars,1);
b(end,1) = 1;


for jj = 1:num_vars
	
	%have to do the matrix inversion for the nullvector
	v = M( (1:num_vars+1)~=jj,:)\b;
	
	detjac_degree = function_degrees(jj)-1;  %this line is correct because the system is homogeneous.
	
	detjac_start = sym(1);
	for ii = 1:detjac_degree
		target_line_coefficients = rand_complex(num_vars+1,1);
		% do a linear move, collect the solutions
		new_points = linear_move(linmove_input, initial_points, random_line_coefficients, target_line_coefficients,'nonsingular_solutions');
		
		% pad the bottom of the new points with the new v values, where v
		% is the nullvector found by matrix inversion
		for qq = 1:size(new_points,2)
			concatenated_point = cell(2*num_vars,1);
			concatenated_point(1:num_vars) = new_points{1,qq};
			for ww = 1:num_vars
				concatenated_point{num_vars+ww} = sprintf('%1.16f %1.16f',real(v(ww)), imag(v(ww)));
			end
			W.add_point(concatenated_point);
		end
		

		
		detjac_input.declare_and_define(sprintf('slice_constant_%i_%i',jj,ii), target_line_coefficients(1), 'constant');

		for kk = 1:num_vars
			detjac_input.declare_and_define(sprintf('LIN_COEFF_%i_%i_%i',kk,jj,ii), target_line_coefficients(kk+1), 'constant');
		end

		line = sym(sprintf('slice_constant_%i_%i',jj,ii));
		for kk=1:length(system_variables)
			line = line + sym(sprintf('LIN_COEFF_%i_%i_%i',kk,jj,ii)) * sym(system_variables{kk});
		end
		detjac_input.declare_and_define(sprintf('line_%i_%i',jj,ii),line,'subfunction');
		detjac_start = detjac_start*sym(sprintf('line_%i_%i',jj,ii));
		
	end
	
	line = sym(0);
	for kk=1:length(system_variables)
		line = line + sym(sprintf('M_%i_%i',jj,kk)) * sym( sprintf('NULLVECTOR_%i',kk) );
	end
	detjac_input.declare_and_define(sprintf('NULLline_%i',jj),line,'subfunction');
	detjac_start = detjac_start*sym(sprintf('NULLline_%i',jj));
	
	
	detjac_input.declare_and_define(sprintf('detjac_start_%i',jj),detjac_start,'subfunction');
	detjac_input.declare_and_define(sprintf('detjac_%i',jj),sprintf('s*detjac_start_%i + (1-s)*detjac_target_%i',jj,jj),'function');
end


line = sym(-1);
jj = num_vars+1;
for kk=1:length(system_variables)
	line = line + sym(sprintf('M_%i_%i',jj,kk)) * sym( sprintf('NULLVECTOR_%i',kk) );
end
detjac_input.declare_and_define('NULLVECTOR_PATCH',line,'function');
	



display('done regenerating start solutions.');


generic_start_points = W.points;
saveme.previous_b_input = b_input_copyme;
saveme.previous_generic_start_points = generic_start_points;
saveme.previous_detjac_input = detjac_input; %#ok 
save('successful_run.mat','-struct','saveme');




cd ..

end













function [detjac_input] = compute_detjac_function(detjac_input)

fprintf('preparing critical points input file\n');


%
% this adds the following fields to the detjac_input struct
%
%  constant.DETJAC_COEFF_varname
%  subfunction.(sprintf('DIFF_%s_%s'))
%
% detjac is a vector of functions, one for each of the original
% functions, plus one for the random combination of the original variables.
%

system_variables = detjac_input.variable_group;
num_vars = length(system_variables);


generic_coefficients = rand_complex(1,num_vars);





vars = sym(zeros(1,num_vars));
for ii = 1:num_vars
	eval(sprintf('syms %s',system_variables{ii}));  %make the variable be a symbol in matlab memory
	vars(ii) = sym(system_variables{ii}); 
	%vars is used in the jacobian call, so that only those derivatives are computed
	detjac_input.declare_and_define(sprintf('DETJAC_COEFF_%s',system_variables{ii}),generic_coefficients(ii),'constant');
end


constant_names = detjac_input.constant(:,1);
for ii = 1:length(constant_names)
	eval(sprintf('syms %s',constant_names{ii}));
end



% the following block of code creates the subfunctions as functions of the
% variables.
%
if ~isempty(detjac_input.subfunction)
	subfunc_names = detjac_input.subfunction(:,1);
else
	subfunc_names = {};
end
num_subfuncs = length(subfunc_names);

if num_subfuncs > 0
	[var_deps,subfunc_deps] = dependency_graph(detjac_input); % determine which subfunctions depend on what, computed by regular expressions
end


absolute_var_deps = zeros(num_subfuncs, num_vars);
for ii= 1:num_subfuncs
	for jj = 1:num_vars
		absolute_var_deps(ii,jj) = depends_on_var(ii,jj, var_deps, subfunc_deps);
	end
end



orig_subfunc_args = cell(num_subfuncs,2);

fprintf('\tsetting up subfunctions in memory\n')

for ii = 1:num_subfuncs
	
	curr_subfunc_args = '';
	curr_subfunc_args_for_regexp = '';
	for jj = 1:num_vars
		if absolute_var_deps(ii,jj)
			curr_subfunc_args = sprintf('%s%s, ',curr_subfunc_args,system_variables{jj});
			curr_subfunc_args_for_regexp = sprintf('%s%s,\\s*',curr_subfunc_args_for_regexp,system_variables{jj});
		end
	end
	curr_subfunc_args = curr_subfunc_args(1:end-2);
	curr_subfunc_args_for_regexp = [curr_subfunc_args_for_regexp(1:end-4) '\s*'];
	
	
	orig_subfunc_args{ii,1} = curr_subfunc_args; %store it.
	orig_subfunc_args{ii,2} = curr_subfunc_args_for_regexp; %store it.
	
	
	currstr = sprintf('syms %s(%s)',detjac_input.subfunction{ii,1},curr_subfunc_args);
	eval(currstr);
end



%this block creates the functions, preserving subfunctions as subfunctions,
%without substitution.

fprintf('\tmaking functions in memory\n');
function_names = detjac_input.functions(:,1);
f = sym(zeros(length(function_names)+1,1));
for ii = 1:length(function_names)
	curr_func = detjac_input.functions{ii,2};
	for jj = 1:num_subfuncs
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{jj});
		newname = sprintf('%s(%s)',subfunc_names{jj},orig_subfunc_args{jj,1});
		newpattern = sprintf('$1%s$2',newname);
		curr_func = regexprep(curr_func,oldpattern,newpattern); 
	end
	f(ii) = eval(curr_func);
end


f(end) = sym(0);
for ii = 1:num_vars
	f(end) = f(end)+sym(sprintf('DETJAC_COEFF_%s',system_variables{ii}))*sym(system_variables{ii});
end



%computes the determinant of the jabobian of the functions with respect to
%the variables.  this respects the chain rule, and produces diff()
%statements in the resulting expression.

fprintf('\tcomputing jacobian matrix\n')


J = jacobian(f,vars);
save('J.mat','J');

v = sym(zeros(num_vars,1));
for ii = 1:num_vars
	eval(sprintf('syms NULLVECTOR_%i',ii));
	v(ii) = eval(sprintf('NULLVECTOR_%i',ii));
end


detjac = J*v;  %this variable is used in some constructed statements below. 




% initialize to empty
new_subfunc_names = {};

%now we substitute away the chain rule diff() statements, by the names of
%new subfunctions, which will be the appropriate partial derivatives.




is_zero_derivative = zeros(length(subfunc_names),num_vars);

fprintf('\tcomputing partial derivatives of subfunctions\n');

for ii = 1:length(subfunc_names)
	for jj = 1:num_vars
		
		if absolute_var_deps(ii,jj)~=0
			
			curr_subfunc = detjac_input.symbol_value(subfunc_names{ii});
			
			for kk = 1:num_subfuncs
				oldpattern = sprintf('(\\W|^)%s(\\W|$)',subfunc_names{kk});
				newname = sprintf('%s(%s)',subfunc_names{kk},orig_subfunc_args{kk,1});
				newpattern = sprintf('$1%s$2',newname);
				curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern); 
			end
			
			
			current_derivative = diff(eval(curr_subfunc),system_variables{jj});
			if current_derivative==0
				is_zero_derivative(ii,jj) = 1;
			end
			
			detjac_input.declare_and_define((sprintf('DIFF_%s_%s',subfunc_names{ii},system_variables{jj})),current_derivative,'subfunction');
			new_subfunc_names{end+1} = sprintf('DIFF_%s_%s',subfunc_names{ii},system_variables{jj});
		else
			is_zero_derivative(ii,jj) = 1;
		end
		
		
	end
end



fprintf('\tdoing partial derivative substitutions for detjac\n')

for ii = 1:length(detjac)
	curr_detjac = char(detjac(ii));
	for jj = 1:length(subfunc_names)
		base = sprintf('diff\\(%s\\(%s\\)',subfunc_names{jj},orig_subfunc_args{jj,2});
		for kk = 1:num_vars

			oldname = sprintf('%s,\\s*%s\\)',base,system_variables{kk});
			oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
			newname = sprintf('DIFF_%s_%s',subfunc_names{jj},system_variables{kk});
			newpattern = sprintf('$1%s$2',newname);
			curr_detjac = regexprep(curr_detjac,oldpattern,newpattern);
		end
	end
	detjac(ii) = curr_detjac;
end





fprintf('\tdoing subfuncion substitutions for detjac\n')

%now we substitute away the subfunction f(...) statements, by the names of
%the original subfunctions.


for ii = 1:length(detjac)
	curr_detjac = char(detjac(ii));
	for jj = 1:length(subfunc_names)
		oldname = sprintf('%s\\(%s\\)',subfunc_names{jj},orig_subfunc_args{jj,2});
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
		newname = sprintf('%s',subfunc_names{jj});
		newpattern = sprintf('$1%s$2',newname);
		curr_detjac = regexprep(curr_detjac,oldpattern,newpattern);
	end
	detjac(ii) = curr_detjac;
end





% finally, we need to substitute away any remaining partial derivative
% statements in the subfunctions.
orig_subfunc_names = subfunc_names;

fprintf('\tdoing substitutions for new subfunctions\n');

for ii = length(new_subfunc_names):-1:1
	ind = detjac_input.symbol_index(new_subfunc_names{ii},'subfunction');
	curr_subfunc = char(detjac_input.subfunction{ind,2});
	for jj = length(orig_subfunc_names):-1:1
		for kk = 1:num_vars
			oldname = sprintf('diff\\(%s\\(%s\\),\\s*%s)',orig_subfunc_names{jj},orig_subfunc_args{jj,2},system_variables{kk});
			oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
			
			newname = sprintf('DIFF_%s_%s',orig_subfunc_names{jj},system_variables{kk});
			newpattern = sprintf('$1%s$2',newname);
			curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
		end
		
		oldname = sprintf('%s\\(%s\\)',orig_subfunc_names{jj},orig_subfunc_args{jj,2});
		oldpattern = sprintf('(\\W|^)%s(\\W|$)',oldname);
		
		newname = orig_subfunc_names{jj};
		newpattern = sprintf('$1%s$2',newname);
		curr_subfunc = regexprep(curr_subfunc,oldpattern,newpattern);
	end
	detjac_input.subfunction{ind,2} = curr_subfunc;
end


for ii = 1:num_vars
	detjac_input.declare_symbols({sprintf('NULLVECTOR_%i',ii)},'variable_group');
	detjac_input.declare_and_define(sprintf('detjac_target_%i',ii),detjac(ii),'subfunction');
end


end %compute_detjac_function








function degrees = get_function_degrees(b_input_copyme)

b_input = copy(b_input_copyme);

dirname = 'degree_computation';
if exist(dirname,'dir')==0
	mkdir(dirname)
end

cd(dirname)


b_input.config.tracktype = '-4';
b_input.config.userhomotopy = '0';
b_input.config.deletetempfiles = '0';

system_variables = b_input.variable_group;
num_vars = length(system_variables);


write_generic_solns(rand_complex(num_vars,1),'start');


write_bertini_input_file(b_input.variable_group, b_input.functions,'filename','input_deg','options',b_input.config,'constants',b_input.constant,'subfunctions',b_input.subfunction);


bertini('filename','input_deg','startname','start','stifle','failcheck',false);

fid = fopen('deg.out','r');

degrees = fscanf(fid,'%i',[num_vars-1 1]);

fclose(fid);

cd ..


end








