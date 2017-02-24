for ii = 1:num_vars
	subfunctions.([system_variables{ii} '_lin']) = sprintf('(1-s)*%s_here + s*%s_rand',system_variables{ii},system_variables{ii});
end
subfunctions.('slice_constant') = '(1-s)*slice_constant_here + s*slice_constant_rand';
