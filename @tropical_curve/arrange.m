function [] = arrange(this, plot_indices)
%% ARRANGE
%
% plot a real tropical curve, arranging the rays according to which
% quadrants had the intersection points, and which sides of the axes (plus
% / minus) had real paths passing through the intersection points.
%
% for best results, dehomogenize the curve first.
%
%  example: arrange(T);
%
% arrange(T,[2 4 5]); arrange onto the 2,4,5 variables.
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



if nargin<2
	if this.num_vars == 2
		plot_indices = [1 2];
	elseif this.num_vars==3
		plot_indices = [1 2 3];
	else
		display('using default variable indices, [1 2 3]');
		plot_indices = [1 2 3];
	end
end

num_plotting_vars = length(plot_indices);

fontsize = 18;
width = 3;


if ~this.isreal()
	error('arrange() is intended for use on real tropical curves.');
end



fig = gcf; %handle to the figure into which we will render.

%the number of rays total in the tropicalization
num_rays = size(this.rays,2);

ray_origin = max(-this.rays,[],2);
ray_origin = max(ray_origin,ones(size(ray_origin)));

% this is for the axes box.
box = max(abs(this.rays),[],2)+1;

%preallocate some colors
colors = parula(num_rays+2);
%discard the first and last.  i think they're too much.
colors = colors(2:end-1,:);

axis_bounds = zeros(num_plotting_vars,2);
for ii = 1:num_rays
	
	%unpack a little
	path_indices = this.aux_data.corresponding_paths{ii};
	
	ray = this.rays(:,ii);
	
	
	%each ray has a multiplicity, and there are that many paths leading to
	%the center point...  we have about three points on each path.  the
	%center, and the pos and neg slices, although a path doesn't have to
	%have points on both slices, just at least one.
	for jj = 1:length(path_indices)
		path_index = path_indices(jj);
		
		p = this.aux_data.paths(path_index);
		v = p.path_variable;
		CP = real(p.connecting_points);
		
		center = real(p.center_point);
		
		for kk = 1:size(CP,2)
			
			connecting_point = CP(:,kk);
			
			
% 			vec = sign(center).*this.rays(:,ii);
			% 			
			ori = ray_origin .* sign(CP(:,kk));
			ray = sign(CP(:,kk)).*this.rays(:,ii);
			switch num_plotting_vars
				case 2
					h = quiver(ori(plot_indices(1)),ori(plot_indices(2)),...
								ray(plot_indices(1)),ray(plot_indices(2)),'MaxHeadSize',1/norm([ray(plot_indices(1)) ray(plot_indices(2))]));
				case 3
					h = quiver3(ori(plot_indices(1)),ori(plot_indices(2)),ori(plot_indices(3)),...
						ray(plot_indices(1)),ray(plot_indices(2)),ray(plot_indices(3)),...
						'MaxHeadSize',1/norm([ray(plot_indices(1)) ray(plot_indices(2))  ray(plot_indices(3))]));
					
					zlabel(sprintf('$%s$',this.variables{plot_indices(3)}),'Interpreter','latex')
				otherwise
					error('unable to arrange without 2 or 3 plotting indices')
			end
			
			xlabel(sprintf('$%s$',this.variables{plot_indices(1)}),'Interpreter','latex')
			ylabel(sprintf('$%s$',this.variables{plot_indices(2)}),'Interpreter','latex')
			axis_bounds = [min(axis_bounds(:,1),ori(plot_indices)+ray(plot_indices)) max(axis_bounds(:,2),ori(plot_indices)+ray(plot_indices))];
			
			set(h,'Color',colors(ii,:),'LineWidth',width);
			hold on
		end
	end
	
	
	

end


if ~isempty(which('labelText'))
	labelText('fontsize',fontsize,'FontName','Times New Roman');
else
	display('note: automatic text resizing available if command ''labelText'' is installed');
end

axis_bounds = [axis_bounds(:,1)-1 axis_bounds(:,2)+1];

if num_plotting_vars
	axis([axis_bounds(1,1) axis_bounds(1,2) axis_bounds(2,1) axis_bounds(2,2)]);
else
	axis([axis_bounds(1,1) axis_bounds(1,2) axis_bounds(2,1) axis_bounds(2,2) axis_bounds(3,1) axis_bounds(3,2)]);
end


% set(gca,'DataAspectRatio',[1 1 1])
hold off
end


