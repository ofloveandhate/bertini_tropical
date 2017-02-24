function plot_cauchy_walks(this, varargin)
error('this function is broken because the place the paths are stored has changed.');
if this.num_vars > 2
	error('plotting the cauchy walks not supported for #vars>3')
end


walked_variables = fieldnames(this.aux_data.paths);
for ii = 1:size(walked_variables)
	
	for jj = 1:length(this.aux_data.paths.(walked_variables{ii}).walked_paths)
		curr_path = this.aux_data.paths.(walked_variables{ii}).walked_paths(jj);
		curr_path.plot();
		hold on
	end
	
end


end
