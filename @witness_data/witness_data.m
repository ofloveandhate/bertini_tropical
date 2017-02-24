% classdef witness_data < handle
%
%  provides a witness_data type for Bertini numerical irreducible
%  decompositions.
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


classdef witness_data < handle
	properties (SetAccess = public)
		num_vars;
		nonempty_codimensions;
		witness_sets;
		
		dehomogenized = false;
	end %re: properties
	
	methods
	
		function obj = witness_data(varargin)
			
			dehom = 0;
			for ii = 1:length(varargin)
				curr_opt = varargin{ii};
				
				switch curr_opt
					case 'parse'
						obj = parse(obj,'witness_data');
						
					case 'filename'
						obj = parse(obj, varargin{ii+1});
						ii = ii+1; %i know, the index is changed in the loop.
						
					case 'dehom'
						dehom = 1;
					otherwise
						error('bad option %s',curr_opt);
				end
			end
			
			if dehom
				obj = dehomogenize(obj);
			end
		end
	
		
		
		
		
		function data = dehomogenize(data)
			if data.dehomogenized
				error('trying to dehomogenize already-dehomogenized witness data.');
			end
			
			for ii = 1:data.num_sets()
				dehomogenize(data.witness_sets(ii)); % pass the call on.
			end
			data.num_vars = data.num_vars-1;
		end
		
		
		function this = reset(this)
			this.num_vars = 0;
			this.nonempty_codimensions = [];
			this.witness_sets = [];
			this.dehomogenized = false;
		end
		
		
		
		function n = num_dims(this)
			n = length(dimensions(this));
		end
		
		
		function d = dimensions(this)
			encountered_dims = zeros(1,this.num_sets());
			for ii = 1:this.num_sets()
				encountered_dims(ii) = this.witness_sets(ii).dimension();
			end
			
			d = unique(encountered_dims);
		end
		
		
		
		function n = num_sets(this)
			n = length(this.witness_sets);
		end
		
		
		function w = choose_sets(this)
			
			
			if this.num_dims ==0
				error('empty witness_data.');
			elseif this.num_dims==1
				fprintf('the witness_data is equidimensional of dimension %i\n\n',this.dimensions());
				d = this.dimensions();  %only one, so just set.
			else
				display('here are the available dimensions:')
				display(this.dimensions()); %print them to the screen
				dims = this.dimensions(); % grab the dimensions
				d = -1;  %initialize to impossible value
				while (isempty(find(dims==d,1)))
					d = input('select a dimension: ');
				end
			end
			
			%get the user's choice of components
			selections = unidimensionsal_select(this, d);
			
			w = construct_agglomeration(this,selections);
		end
		
		
		
		function w = construct_agglomeration(this,indices)
			
			w = witness_set(); %make an empty witness set.
			
			
			for ii = 1:length(indices)
				curr_set = indices(ii);
				
				w.add_point(this.witness_sets(curr_set).points);
				
				w.component_number = [w.component_number this.witness_sets(curr_set).component_number];
				if ii==1
					w.add_linear(this.witness_sets(curr_set).linears);
					w.add_patch(this.witness_sets(curr_set).patches);
					w.set_randomization_matrix(this.witness_sets(curr_set).randomization_matrix);
					w.codimension = this.witness_sets(curr_set).codimension;
					w.singular = this.witness_sets(curr_set).singular;
					w.dehomogenized = this.witness_sets(curr_set).dehomogenized;
					
				else
					if w.singular ~= this.witness_sets(curr_set).singular;
						error('mixing singular and nonsingular components is not permitted');
					end
				end
			end %re: for ii
			
		end %re: construct_agglomeration
		
		
		
		
		function selections = unidimensionsal_select(this, dim)
			
			
			matches = 0;
			for ii = 1:this.num_sets()
				if this.witness_sets(ii).dimension==dim
					matches = matches + 1;
				end
			end
			
			
			if matches==0
				error('there are no components of the requested dimension (%i).',dim);
			elseif matches==1
				tmpselections = 1;
			else
				fprintf('\nhere are the available components:\n\n');

				matches = 0;
				for ii = 1:this.num_sets()
					if this.witness_sets(ii).dimension==dim
						matches = matches + 1;
						fprintf('%i:  degree %i, singular %i\n',matches, this.witness_sets(ii).degree, this.witness_sets(ii).singular);
					end
				end
				fprintf('\n');


				ncomponents = input('how many components would you like to select? ');
				
				if ncomponents~=matches
					tmpselections = [];

					while length(tmpselections) < ncomponents
						c = input('component selection: ');
						tmpselections = unique([tmpselections c]);
					end				
				else
					tmpselections = 1:matches;
				end
			end
			
			selections = [];
			component_count = 0;
			for ii = 1:this.num_sets()
				if this.witness_sets(ii).dimension == dim
					component_count = component_count+1;
					if ~isempty(find(tmpselections==component_count,1))
						 selections(end+1) = ii; %thanks, but i don't want to preallocate.
					end
				end
			end
		end % unidimensional_select
		
		%FUNCTIONS DEFINED EXTERNALLY
		data = parse(data, filename)
	
		
		
		
		
	end %re: methods
end %re: classdef
