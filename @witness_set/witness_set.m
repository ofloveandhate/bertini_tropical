% classdef witness_set  < handle
%
%  Provides a witness_set class, for interacting with witness_data files
%  from Bertini 1.
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


classdef witness_set  < handle
	properties
		points = [];
		linears = [];
		patches;
		
		randomization_matrix = [];
		dehomogenized = false;
		
		component_number;
		codimension;
		
		singular;
		
	end
	
	
	methods
		
		
		function obj = witness_set()
			obj.dehomogenized = false;
		end
		
		function this = add_point(this, new_point)
			this.points = [this.points new_point];
		end
		
		function n = num_points(this)
			n = size(this.points,2);
		end
		
		
		
		function d = degree(this)
			d = size(this.points,2);
		end
		
		
		function n = num_variables(this)
			n = size(this.points,1);
		end
		
		function this = add_linear(this, new_linear)
			this.linears = [this.linears new_linear];
		end
		
		function n = num_linears(this)
			n = size(this.linears,2);
		end
		
		function this = add_patch(this, new_patch)
			this.patches = [this.patches new_patch];
		end
		
		function this = set_randomization_matrix(this, R)
			if ~ismatrix(R)
				error('trying to set non-matrix object as randomization matrix for witness_set')
			end
			
			this.randomization_matrix = R;
			
		end
		
		function t = is_consistent(this)
			if this.dehomogenized
				t = size(this.linears,2) == (size(this.points,1) - this.codimension);
			else
				t = size(this.linears,2) == (size(this.points,1) - this.codimension - 1);
			end
		end

		
		function d = dimension(this)
			if ~this.is_consistent()
				error('inconsistent witness set, the number of linears does not jive with the number other things');
			end
			
			d = this.num_linears();
			
		end
		
		function this = dehomogenize(this)
			if this.dehomogenized
				error('trying to dehomogenize an already-dehomogenized witness set.');
			end
			
			this.points = dehomogenize(this.points,1,1); 
			%use the first dimension (rows) for dehomogenization.  use the first coordinate as the (de)homogenizing coordinate
			
			% change the flag over
			this.dehomogenized = true;
			
		end
		
		
		
	end
	
end %re: classdef witness_set
