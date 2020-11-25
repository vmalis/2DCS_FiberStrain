function [minIntrVec,stat,actpctg] = genSampling(pdf,iter,tol)

%  [mask,stat,N] = genSampling(pdf,iter,tol)
%
%  UC San Diego / March 2019 / Vadim Malis

%h = waitbar(0);

pdf(find(pdf>1)) = 1;
K = sum(pdf(:));

minIntr = 1e99;
minIntrVec = zeros(size(pdf));

for n=1:iter
	tmp = zeros(size(pdf));
	while abs(sum(tmp(:)) - K) > tol
		tmp = rand(size(pdf))<pdf;
	end
	
	TMP = ifft2(tmp./pdf);
	if max(abs(TMP(2:end))) < minIntr
		minIntr = max(abs(TMP(2:end)));
		minIntrVec = tmp;
	end
	stat(n) = max(abs(TMP(2:end)));
	%waitbar(n/iter,h);
end

actpctg = sum(minIntrVec(:))/prod(size(minIntrVec));

%close(h);


