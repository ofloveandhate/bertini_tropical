% class for holding paths in complex time around the origin.
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

classdef cauchy_path < handle
	properties
		cycle_number;  %how many times you have to go around the origin.
		path;          %the computed path around the origin
		start_point;   %the start point for the path.  this seems redundant.
		path_variable; %which variable was used to go around the origin.
		radius;
		center_point;
		connecting_points;
		
		aux_data; % any miscellaneous data computed about the path.
	end
	
	
	
	methods
		
		
		function obj = add_connecting_point(obj, new_point)
			if isempty(obj.connecting_points)
				obj.connecting_points = new_point;
			else
				obj.connecting_points = [obj.connecting_points new_point];
			end
		end
		
		function obj = cauchy_path()
			obj.cycle_number = 0;
			obj.path = [];
			obj.start_point = [];
		end
		
		function obj = set_start_point(obj, new_point)
			obj.path(:,1) = new_point;
			obj.start_point = new_point;
		end
		
		
		function obj = add_point_to_path(obj, new_point)
			obj.path(:,end+1) = new_point;
		end
		
		function obj = set_center(obj, new_center)
			obj.center_point = new_center;
		end
		
		function obj = set_radius(obj,new_radius)
			obj.radius = new_radius;
		end
		
		function disp(obj)
			fprintf('cycle number: %i\n\n',obj.cycle_number)
			disp('the followed path:');
			disp(obj.path);
			disp('center point:')
			disp(obj.center_point)
			disp('connecting points:');
			disp(obj.connecting_points);
		end
		
		
		function this = dehomogenize(this, homvar_index)
			
			nvars = this.num_vars();
			this.path = this.path(1:nvars~=homvar_index,:)./repmat(this.path(homvar_index,:),[nvars-1 1]);
			this.center_point = this.center_point(1:nvars~=homvar_index)/this.center_point(homvar_index);
			
			this.connecting_points = this.connecting_points(1:nvars~=homvar_index,:)./repmat(this.connecting_points(homvar_index,:),[size(this.connecting_points,1)-1 1]);
		end
		
		
		function this = rescale(this)
			this.path = rescale(this.path);
			this.center_point = rescale(this.center_point);
			this.connecting_points = rescale(this.connecting_points);
		end
		
		function this = unscale(this, unscaling_function)
			
			this.path = unscaling_function(this.path);
			this.center_point = unscaling_function(this.center_point);
			this.connecting_points = unscaling_function(this.connecting_points);
		end
		
		function n = num_vars(this)
			n = size(this.path,1);
		end
		
		function h = plot(this,varargin)
			
			if this.num_vars()~=2
				error('plotting cauchy paths not supported for more than 2 variables')
			end
			
			scaled_path = this.path;
% 			scaled_path = real(scaled_path) + 10*imag(scaled_path);
			h = [];
			
			
			x = real(scaled_path(1,:));
			y = real(scaled_path(2,:)) ;
			
			z = imag(scaled_path(1,:));
			
			
			
			h(1) = plot3(x, y, z);
% 			hold on
% 			
% 			z = imag(scaled_path(2,:));
% 			
% 			h(2) = plot3(x, y, z);
			
			h(end+1) = plot3(real(this.center_point(1)), real(this.center_point(2)), 0,'Color','r','Marker','x');
			
		end
	end
	
	
end
