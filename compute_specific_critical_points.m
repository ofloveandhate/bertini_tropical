% critpts = compute_specific_critical_points(b_input_copyme,var_index,generic_critical_points)
%
% move generic coordinates for the jacobian to particular ones for the
% specific variable being considered as time.
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


function critpts = compute_specific_critical_points(b_input_copyme,var_index,generic_critical_points)


b_input = copy(b_input_copyme);


system_variables = b_input.variable_group;
num_vars = length(system_variables);





dirname = sprintf('critpts_%s',system_variables{var_index});
if exist(dirname,'dir')==0
	mkdir(dirname)
end
cd(dirname)


if exist('successful_run.mat','file')
	load('successful_run.mat');
	
	if issame_system(previous_b_input,b_input)
		display('using previously computed critpts');
		critpts = previous_critpts;
		cd ..
		return
	end
end


fprintf('computing critical points for t = %s\n',system_variables{var_index});




for ii = 1:num_vars/2
	if ii==var_index
		b_input.define_symbol(sprintf('DETJAC_COEFF_%s',system_variables{ii}),'1');
	else
		b_input.define_symbol(sprintf('DETJAC_COEFF_%s',system_variables{ii}),'0');
	end
end

write_generic_solns(generic_critical_points,'start_critpts');


% the degreebound number probably changed.
if isfield(b_input.config,'degreebound')
	b_input.config = rmfield(b_input.config,'degreebound');
end


b_input.config.tracktype = '0';
b_input.config.userhomotopy = '2';

write_generic_solns(generic_critical_points,'start_critpts');

evalme = 'write_bertini_input_file(system_variables, b_input.functions,''filename'',''input_specific_critpts'',''options'',b_input.config';
evalme = sprintf('%s,''subfunctions'',b_input.subfunction',evalme);
evalme = sprintf('%s,''constants'',b_input.constant',evalme);
evalme = sprintf('%s);',evalme);
eval(evalme);


bertini('filename','input_specific_critpts','startname','start_critpts','stifle','mpi',min(4,size(generic_critical_points,2)));

critpts = get_generic_solns('finite_solutions',num_vars);

critpts = critpts(1:end/2,:); % trim off the synthetic variables.


saveme.previous_b_input = b_input_copyme;
saveme.previous_critpts = critpts; %#ok 
save('successful_run.mat','-struct','saveme');

cd ..

end
