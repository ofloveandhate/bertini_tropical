% status = bertini(varargin)
%
% a wrapper for calling command-line software
% bertini from Matlab.  returns the returned value of the system call to
% Bertini.
%
% valid options to bertini() are:
% 
% * 'mpi', with a second trailing argument being the number of processes to
% use.  This of course requires MPI to be installed on your machine.
% * 'stifle' -- Redirects the screen output to a named file, automatically
% named as the time of the run.
% * 'nostifle' -- Turns off redirecting of screen output.
% * 'startname' -- The name of the start file to use, if needed.
% * 'filename' -- The name of the input file.
% * 'showcommand' -- Display the system command before calling it.
% * 'failcheck' -- Check the file 'failed_paths' for failed paths.
% * 'rerunonfail' -- If there are failed paths, re-run with the same
% settings.  This can be useful, because there are random numbers behind
% the scene in Bertini which can influence the outcome of the run, for
% better or for worse.
%
% copyright 2015, 2016, 2017 Dani Brake
% University of Notre Dame
% Applied and Computational Mathematics and Statistics
% University of Wisconsin, Eau Claire
% Mathematics
% danielthebrake@gmail.com brakeda@uwec.edu
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


function status = bertini(varargin)


rerun_on_fail = 0;
failcheck = false;
inputfile = 'input';
startfile = '';
stifletext = '';
mpitext = '';
showcommand = false;

ii = 1;
while ii<=length(varargin)
	switch varargin{ii}
		case 'mpi'
			numprocs = varargin{ii+1};
			if ~isint(numprocs)
				error('argument trailing ''mpi'' -- numprocs -- must be a positive integer');
			end
			if numprocs <= 0
				error('argument trailing ''mpi'' -- numprocs -- must be a positive integer');
			end
			
			if and(numprocs > 1, system('which mpiexec > /dev/null')==0)
				mpitext = sprintf('mpiexec -n %i',numprocs);
			end
			ii = ii+1;
			
		case 'stifle'
			stifletext = '> /dev/null';
% 			stifletext = sprintf('> run_results_%s',datestr(now,'yymmddHHMMSSfff'));
		
		case 'nostifle'
			stifletext = '';
			
		case 'startname'
			startfile = varargin{ii+1};
			ii = ii+1;
			
		case 'filename'
			
			inputfile = varargin{ii+1};
			ii = ii+1;
			
		case 'showcommand'
			showcommand = true;
		case 'failcheck'
			failcheck= varargin{ii+1};
			ii = ii+1;
			
		case 'rerunonfail'
			failcheck = true;
			rerun_on_fail= varargin{ii+1};
			ii = ii+1;
			
		otherwise
			error('bad option ''%s'' to bertini ()',varargin{ii})
	end
	
	ii = ii+1;
end


command = sprintf('%s bertini %s %s %s',mpitext, inputfile,startfile,stifletext);

if showcommand
	display(command);
end

status = system(command);

if failcheck
	
	
	[had_failures, failure_messages] = parse_failed_paths();
	
	if had_failures
		warning('had failed paths in this call of bertini');
		
		
		if rerun_on_fail>0
			warning('attempting to run again with same settings...');
			
			successful_resolution = false;
			for ii = 1:rerun_on_fail
				status = system(command);

				[had_failures, failure_messages] = parse_failed_paths();

				if ~had_failures
					successful_resolution = true;
					break
				end;
			end

			if successful_resolution
				display('resolved');
			else
				warning('continued failure.  adjust bertini settings.');
			end
		end
		
		

	end
	
	
	
end




if status~=0
	error('bertini was not successful');
end
end


function [had_failures, failure_messages] = parse_failed_paths()

had_failures = false;
failure_messages = {};


fid = fopen('failed_paths','r');
topline = fgetl(fid);

% feof(fid);


fclose(fid);

if ~isempty(topline)
	had_failures = true;
end

end


