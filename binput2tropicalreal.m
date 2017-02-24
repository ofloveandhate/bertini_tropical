%  [bertini_options, system_variables, curve_system, subfunctions, is_homogenized, homogenizing_variable] = binput2tropicalreal(b_input)
%
%  Take a bertini input file (matlab representation from this package), and
%  extract data from it for the tropical decomposition.
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


function [bertini_options, system_variables, curve_system, subfunctions, is_homogenized, homogenizing_variable] = binput2tropicalreal(b_input)


% if length(b_input.variable_group)~=1
% 	error('sorry, tropical real is only enabled for exactly one variable group');
% end


system_variables = b_input.variable_group;

curve_system = b_input.functions;
bertini_options = b_input.config;
subfunctions = b_input.subfunction;

if isfield(bertini_options,'tracktype')
	bertini_options = rmfield(bertini_options,'tracktype');
end


found_h = false;
for ii = 1:length(system_variables)
	if system_variables{ii}=='h'
		found_h = true;
	end
end


if ~found_h
	is_homogenized = 0;
	homogenizing_variable = '';
% 	error('tropical real expects there to be a variable named ''h'', and for it to act as the homogenizing variable.');
else
	homogenizing_variable = 'h';
	is_homogenized = 1;
end






end
