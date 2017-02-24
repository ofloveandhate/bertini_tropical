function g = gcd_all(A)
% compute the gcd of all inputs.  does so recursively.  if any are 0

A = abs(A(A~=0));

g = A(1);
for k = 2:length(A)
	 g = gcd(g,A(k));
	 if g==1
		 return;
	 end
end


end
