% connecting_point = connect_slice_to_intersection(b_input_copyme, start_point, start_constant, target_constant, variable_index, reality_string)
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

function connecting_point = connect_slice_to_intersection(b_input_copyme, start_point, start_constant, target_constant, variable_index, reality_string)

if ~or(strcmp(reality_string,'real'),strcmp(reality_string,'complex'))
	error('reality_string must be either ''real'' or ''complex''');
end

b_input = copy(b_input_copyme);

system_variables = b_input.variable_group;

curr_varname = system_variables{variable_index};



num_vars = length(start_point);



b_input.config.userhomotopy = '2';


b_input.declare_and_define('slice_constant_rand',start_constant,'constant');
b_input.declare_and_define('slice_constant_here',target_constant,'constant');
b_input.declare_and_define('slice_constant','(1-s)*slice_constant_here + s*slice_constant_rand','subfunction');






line = -sym('slice_constant')+genvarname(curr_varname);

b_input.declare_and_define('moving_line',line,'function');



dirname = sprintf('%s_connect',system_variables{variable_index});
if exist(dirname,'dir')==0
	mkdir(dirname)
end
cd(dirname)



write_generic_solns(start_point,'start_connect');


write_bertini_input_file(system_variables, b_input.functions,'filename','input_connect','options',b_input.config,'constants',b_input.constant,'subfunctions',b_input.subfunction);


bertini('filename','input_connect','startname','start_connect','stifle','mpi',min(4,size(start_point,2)),'rerunonfail',5);%

if strcmp(reality_string,'real')
	connecting_point = get_generic_solns('real_finite_solutions',num_vars);
else
	connecting_point = get_generic_solns('finite_solutions',num_vars);
end


cd ..

end
