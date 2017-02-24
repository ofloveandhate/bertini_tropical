% val = issame_system(sys1, sys2)
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
function val = issame_system(sys1, sys2)

if isequal(sys1.variable_group,sys2.variable_group) &&...
	isequal(sys1.hom_variable_group,sys2.hom_variable_group) &&...
	isequal(sys1.constant,sys2.constant) &&...
	isequal(sys1.parameter,sys2.parameter) &&...
	isequal(sys1.subfunction,sys2.subfunction) &&...
	isequal(sys1.functions,sys2.functions)
	val = true;
else
	val = false;
end

end
