%linecount: determines the number of lines in a text file.  

%returns [] if the file does not exist
%        [int] the number of lines to be read in if it does.

%input: string file

% copyright 2013, 2016 Daniel Brake
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

function linecount = getlinecount(file)

if ~ischar(file)
	error('input file name was not a string');
end
	
	
if isempty(dir(file))
	error('requested file ''%s'' does not exist.',file);
end


fid = fopen(file,'r');
linecount = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  linecount = linecount+1;
end

fclose(fid);
end





