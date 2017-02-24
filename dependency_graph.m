%[var_deps,subfunc_deps] = dependency_graph(b_input)
%
%computes dependency relationship information among subfunctions for a
%bertini input file, both for variables and for other subfunctions.
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

function [var_deps,subfunc_deps] = dependency_graph(b_input)


system_variables = b_input.variable_group;
num_vars = length(system_variables);


subfunction_names = b_input.subfunction(:,1);
num_subfuncs = length(subfunction_names);


var_deps = zeros(num_subfuncs,num_vars);

for ii = 1:num_subfuncs
	curr_subfunc = char(b_input.subfunction{ii,2});
	for jj = 1:num_vars
		if ~isempty(regexp(curr_subfunc,sprintf('(\\W|^)%s(\\W|$)',system_variables{jj}),'once'))
% 			fprintf('%s depends on %s\n',curr_subfunc,system_variables{jj})
			var_deps(ii,jj) = 1;
		end

	end
end

subfunc_deps = zeros(num_subfuncs,num_subfuncs);

for ii = 1:num_subfuncs
	curr_subfunc = char(b_input.subfunction{ii,2});
	for jj = 1:num_subfuncs
		if ~isempty(regexp(curr_subfunc,sprintf('(\\W|^)%s(\\W|$)',subfunction_names{jj}),'once'))
% 			fprintf('%s depends on %s\n',curr_subfunc,subfunction_names{jj})
			subfunc_deps(ii,jj) = 1;
		end

	end
end



end
