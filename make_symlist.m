symlist = system_variables;

for ii = 1:num_vars
	symlist{end+1} = [system_variables{ii} '_rand'];
	symlist{end+1} = [system_variables{ii} '_here'];
end

symlist{end+1} = 'slice_constant';
symlist{end+1} = 'slice_constant_rand';
symlist{end+1} = 'slice_constant_here';



