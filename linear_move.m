% move_points = linear_move(b_input_copyme, initial_points, random_line_coefficients, target_line_coefficients,filetoread)
%
% Move points which all lie on a line, to a target line, using Bertini.
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

function move_points = linear_move(b_input_copyme, initial_points, random_line_coefficients, target_line_coefficients,filetoread)

b_input = copy(b_input_copyme);


if nargin<5
	filetoread = 'real_finite_solutions';
end
% make_vars
% make_syms
% make_subfuncs % sets up _lin = _here + _rand for each variable, plus the slice_constant


dirname = sprintf('linear_move_tmp');

if exist(dirname,'dir')==0
	mkdir(dirname);
end

cd(dirname);




write_generic_solns(initial_points,'start');

b_input.config.userhomotopy = '2';


system_variables = b_input.variable_group;
num_vars = length(system_variables);


b_input.declare_symbols({'slice_constant_rand' 'slice_constant_here'},'constant');

b_input.set_symbol('slice_constant_rand',random_line_coefficients(1));
b_input.set_symbol('slice_constant_here',target_line_coefficients(1));


b_input.declare_and_define('slice_constant','(1-s)*slice_constant_here + s*slice_constant_rand','subfunction');



for ii = 1:num_vars
	b_input.declare_symbols({[system_variables{ii} '_lin']},'subfunction');
	
	b_input.set_symbol([system_variables{ii} '_lin'], sprintf('(1-s)*%s_here + s*%s_rand',system_variables{ii},system_variables{ii}) );

	b_input.declare_symbols({[char(system_variables(ii)) '_here'] [char(system_variables(ii)) '_rand']},'constant');
	
	b_input.set_symbol([system_variables{ii} '_here'], target_line_coefficients(ii+1));
	b_input.set_symbol([system_variables{ii} '_rand'], random_line_coefficients(ii+1));
end



line = sym(matlab.lang.makeValidName('slice_constant'));
for ii=1:length(system_variables)
	line = line + sym(sprintf('%s_lin',system_variables{ii})) * sym(system_variables{ii});
end
b_input.declare_symbols({'line'},'function');
b_input.set_symbol('line',char(line));


write_bertini_input_file(b_input.variable_group, b_input.functions,'filename','input_linear_move','options',b_input.config,'constants',b_input.constant,'subfunctions',b_input.subfunction);

bertini('filename','input_linear_move','startname','start','stifle','mpi',min(4,size(initial_points,2))); %,'showcommand'


move_points = get_generic_solns(filetoread,num_vars,'str');



cd ..


end





