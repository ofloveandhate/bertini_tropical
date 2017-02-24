% [unique_solns,indices] = unique_up_to_threshold(test_solns, threshold)
% 
% sorts an array for unique columns, up to a numeric threshold on 2-norm of
% difference
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

function [unique_solns,indices] = unique_up_to_threshold(test_solns, threshold)


%returns the unique solutions, up to 0-norm threshold
[num_vars,num_solns] = size(test_solns);


unique_solns = zeros(num_vars,num_solns); % preallocate

num_unique = 0; %set to zero

curr_num_solns = size(test_solns,2); % how many test solutions remain

were_unique = false(1,num_solns);  %preallocate

curr_indices = 1:num_solns; %preallocate.

while curr_num_solns>0
	
	
	curr_soln = test_solns(:,1); %take the first test solution to be next unique one.
	
	were_unique(curr_indices(1)) = true;
	
	rest = test_solns(:,2:end); %trim off the ones behind the first to become 'rest', the remainder of them
	curr_indices = curr_indices(2:end);
	
	
	num_unique = num_unique+1; %increment this counter
	unique_solns(:,num_unique) = curr_soln; % copy current solution into the unique
	
	R = repmat(curr_soln,[1 size(rest,2)]);
	bigind = or(abs(rest)>1,abs(R)>1);
	diffs = abs(rest - R);
	diffs(bigind) = relative_difference(rest(bigind),R(bigind));

	near_pts = all(diffs < threshold,1);
	
	curr_indices = curr_indices(~near_pts);
	test_solns = rest(:,~near_pts); %use logical indexing to get rid of any in the remainder that are close
	curr_num_solns = size(test_solns,2); %see how many remain.  if 0, then we're done!
end

unique_solns = unique_solns(:,1:num_unique); %trim the fat

indices = 1:num_solns;
indices = indices(were_unique);

end
