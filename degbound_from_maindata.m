% degreebound = degbound_from_maindata()
%
% extract the degree bound from the 'main_data' file from Bertini.
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


function degreebound = degbound_from_maindata()


filename = 'main_data';

num_lines = getlinecount(filename);



fid = fopen(filename);


line_count = 0;
while line_count < num_lines
	line_count=line_count+1;
	currline = fgetl(fid);
	
	start = strfind(currline,'DegreeBound: ');
	if start
		endloc = strfind(currline,';');
		degreebound = str2double(currline(start+13:endloc-1));
		break
	end

	
end

			

fclose(fid);



end
