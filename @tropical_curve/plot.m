% plot(this,varargin)
%
% Plot a tropical curve.
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


function plot(this,varargin)
renderaxes = false;
rendertext = true;
rendermult = false;

fig = gcf;

fontsize = 18;

vec = this.vectors;
num_vec = size(vec,2);

zer = zeros(1,num_vec);
switch(this.num_vars)
	case 2
		
		if renderaxes
			a(1) = quiver([0 0],[0 0],[1 0],[0 1]);
			hold on
			
			set(a,'Color','k','LineWidth',2)
			
		end
		
		
		h = quiver(zer,zer,vec(1,:),vec(2,:));
		set(gca,'DataAspectRatio',[1 1 1])
		
		if rendertext
			xlabel(this.variables{1});
			ylabel(this.variables{2});
			
			textme = cell(1,num_vec);
			
			for ii = 1:num_vec
				if rendermult
					textme{ii} = sprintf('  < %i %i >, mult=%i',vec(1,ii),vec(2,ii), this.multiplicities(ii));
				else
					textme{ii} = sprintf('  < %i %i >',vec(1,ii),vec(2,ii));
				end
			end
			
			
			placement_factor = 0.85;
			text(placement_factor*vec(1,:), placement_factor*vec(2,:),textme,...
				'HorizontalAlignment','left',...
				'FontSize',fontsize);
			
			
		end
		
		axis_vec = [min(vec(1,:))-1 max(vec(1,:))+1 min(vec(1,:))-1 max(vec(2,:))+1];
	case 3
		
		if renderaxes
			a(1) = quiver3([0 0 0],[0 0 0],[0 0 0],[1 0 0],[0 1 0],[0 0 1]);
			hold on
			set(a,'Color','k','LineWidth',2)
		end
		
		
		h = quiver3(zer,zer,zer,vec(1,:),vec(2,:),vec(3,:));
		set(gca,'DataAspectRatio',[1 1 1])
		
		if rendertext
			xlabel(this.variables{1},'FontSize',fontsize);
			ylabel(this.variables{2},'FontSize',fontsize);
			zlabel(this.variables{3},'FontSize',fontsize);
			
			textme = cell(1,num_vec);
			
			for ii = 1:num_vec
				if rendermult
					textme{ii} = sprintf('< %i %i %i >, mult=%i',vec(1,ii),vec(2,ii),vec(3,ii), this.multiplicities(ii));
				else
					textme{ii} = sprintf('< %i %i %i >',vec(1,ii),vec(2,ii),vec(3,ii));
				end
				
				
			end
			
			placement_factor = 0.85;
			text(placement_factor*vec(1,:), placement_factor*vec(2,:), placement_factor*vec(3,:),textme,...
				'HorizontalAlignment','left',...
				'FontSize',fontsize);
			
		end
		axis_vec = [-1 max(vec(1,:))+1 -1 max(vec(2,:))+1 -1 max(vec(3,:))+1];
		
		
		
		
	otherwise
		
end

set(h,'LineWidth',4);
set(gca,'FontSize',fontsize)
axis(axis_vec);
grid off

if rendertext
	title(char(this.system_name),'interpreter','none');
end

render_into_file([ this.system_name '_real_tropical_curve'] )



hold off
end
