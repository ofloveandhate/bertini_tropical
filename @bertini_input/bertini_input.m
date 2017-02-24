% classdef bertini_input < matlab.mixin.Copyable
%
%  Provides an object-oriented Bertini input file for Matlab.
%
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

classdef bertini_input < matlab.mixin.Copyable
	
	properties
		
		filename;
		md5hash;
		
		
		declared_symbols;
		config;
		
		variable_group;
		hom_variable_group;
		
		constant;
		parameter;
		subfunction;
		functions;
		%
	end
	
	
	
	methods
		
		
		
		
		
		function b_input = bertini_input(filename)
			
			if nargin < 1
				b_input.filename = 'unset';
			else
				b_input.filename = filename;
				b_input = parse(b_input,b_input.filename);
			end
			
		end
		
		
		
		%defined externally
		b_input = parse(b_input, filename);
		
		% a wrapper function allowing to parse from an input file contained
		% in a string 
		function b_input = parse_from_string(b_input, file_as_str)
			
			tmp_filename = 'bertini_input_temp_file_for_parsing_deleteme';
			fid = fopen(tmp_filename,'w');
			fprintf(fid,'%s',file_as_str);
			fclose(fid);
			
			b_input = parse(b_input,tmp_filename);
			
			delete(tmp_filename);
			
		end
		
		
		function b_input = define_symbol(b_input,name,value)
			b_input.set_symbol(name,value);
		end
		
		
		function b_input = set_symbol(b_input, name, value)
			
			type = '';
			num_symbols = size(b_input.declared_symbols,1);
			for ii = 1:num_symbols
				if strcmp(name, b_input.declared_symbols{ii,1})
					type = b_input.declared_symbols{ii,2};
					break;
				end
			end
			
			if isempty(type)
				error('desired_symbol ''%s'' not found.  please declare first', name);
			end
			
			
			index = 0;
			for ii = 1:size(b_input.(type),1)
				if strcmp(name, b_input.(type){ii,1})
					index = ii;
					break;
				end
			end
			if index==0
				error('desired_symbol ''%s'' not found.  please declare first', name);
			end
			
			b_input.(type){index,2} = value;
		end
		
		
		
		%Check whether the input file already has a symbol of a given name
		%declared.  true if it is.
		function v = is_symbol_declared(b_input, name)
			num_symbols = size(b_input.declared_symbols,1);
			for ii = 1:num_symbols
				if strcmp(b_input.declared_symbols{ii,1},name)
					v = true;
					return;
				end
			end
			v = false;
		end
		
		
		function b_input = declare_symbols(b_input, names, decl_type)
			for mm = 1:length(names)
				if is_symbol_declared(b_input, names{mm})
					error('previously declared symbol %s being declared again as type %s',names{mm},decl_type);
				end
			end
			
			
			
			if strcmp(decl_type,'function')
				type = 'functions';
			else
				type = decl_type;
			end
			
			num_symbols = size(b_input.declared_symbols,1);
			for mm = 1:length(names)
				b_input.declared_symbols{num_symbols+1,1} = names{mm};
				b_input.declared_symbols{num_symbols+1,2} = type;
				num_symbols = num_symbols+1;
				b_input.(type){end+1,1} = names{mm};
			end
			
			
		end
		
		
		
		function b_input = declare_and_define(b_input,name,value,type)
			b_input.declare_symbols({name},type);
			b_input.set_symbol(name,value);
		end
		
		function val = symbol_value(b_input, name)
			
			type = b_input.symbol_type(name);
			val = b_input.(type){b_input.symbol_index(name,type),2};
			
		end
		
		
		function b_input = delete_symbol(b_input, name)
			[type,decl_index] = b_input.symbol_type(name);
			index = b_input.symbol_index(name,type);
			
			num_symbols = size(b_input.declared_symbols,1);
			b_input.declared_symbols = b_input.declared_symbols(1:num_symbols~=decl_index,:);
			
			b_input.(type) = b_input.(type)(1:size(b_input.(type),1)~=index,:);
		end
		
		
		function b_input = change_symbol_type(b_input, name, new_type)
			if strcmp(new_type,'function')
				new_type = 'functions';
			end
			
			[old_type,decl_index] = b_input.symbol_type(name);
			old_index = b_input.symbol_index(name,old_type);
			
			%change the declared type
			b_input.declared_symbols{decl_index,2} = new_type;
			
			b_input.(new_type)(end+1,:) = b_input.(old_type)(old_index,:);
			
			% delete from old type's set of symbols
			b_input.(old_type) = b_input.(old_type)(1:size(b_input.(old_type),1)~=old_index,:);
			
		end
		
		
		function write(b_input, filename)
			
			error('this function needs to be implemented');
			
		end
		
		
		
		function S = saveobj(this)
			S.filename = this.filename;
			S.md5hash = this.md5hash;
			S.declared_symbols = this.declared_symbols;
			S.config = this.config;
			S.variable_group = this.variable_group;
			S.hom_variable_group = this.hom_variable_group;
			S.constant = this.constant;
			S.parameter = this.parameter;
			S.subfunction = this.subfunction;
			S.functions = this.functions;
			
		end
		
		
		
		
		function [type,decl_index] = symbol_type(b_input, name)
			type = '';decl_index = 0;
			num_symbols = size(b_input.declared_symbols,1);
			for ii = 1:num_symbols
				if strcmp(name, b_input.declared_symbols{ii,1})
					type = b_input.declared_symbols{ii,2};
					decl_index = ii;
					break;
				end
			end
			if isempty(type)
				error('desired_symbol ''%s'' not found.', name);
			end
		end
		
		
		function index = symbol_index(b_input,name,type)
			index = 0;
			for ii = 1:size(b_input.(type),1)
				if strcmp(name, b_input.(type){ii,1})
					index = ii;
					break;
				end
			end
			if index==0
				error('desired_symbol ''%s'' not found.', name);
			end
		end
	end % re: private methods
	
	
	
		
		
	methods(Static)
		function obj = loadobj(s)
			if isstruct(s)
				newObj = bertini_input;
				
				newObj.filename = s.filename;
				newObj.md5hash = s.md5hash;
				newObj.declared_symbols = s.declared_symbols;
				newObj.config = s.config;
				newObj.variable_group = s.variable_group;
				newObj.hom_variable_group = s.hom_variable_group;
				newObj.constant = s.constant;
				newObj.parameter = s.parameter;
				newObj.subfunction = s.subfunction;
				newObj.functions = s.functions;

				obj = newObj;
				
			else
				obj = s;
			end
		end
	end
	
end
