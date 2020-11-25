function x = CSL1NlCg(x0,param)
% 
% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
% 
% Given the acquisition model y = E*x, and the sparsifying transform W, 
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda1 * ||W*x||_1 + lambda2 * TV(x) 
% 
%  ---------------------------------------
%  UC San Diego / Januarry - March 2019 / Vadim Malis
%

fprintf('\n Non-linear conjugate gradient algorithm')
fprintf('\n ---------------------------------------\n')

% starting point
x=x0;

% line search parameters
maxlsiter = 150 ;
gradToll = 1e-3 ;
param.l1Smooth = 1e-15;	
alpha = 0.01;  
beta = 0.6;
t0 = 1 ; 
k = 0;

% compute g0  = grad(f(x))
g0 = grad(x,param);
dx = -g0;

% iterations
while(1)

    
    if param.nite_fft == 0
        param.nite = param.nite_pca;
    else
        param.nite = param.nite_fft;
    end
        
    
    % backtracking line-search
	f0 = objective(x,dx,0,param);
	t = t0;
    f1 = objective(x,dx,t,param);
	lsiter = 0;
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		f1 = objective(x,dx,t,param);
	end

	if lsiter == maxlsiter
		disp('Error - line search ...');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2, t0 = t0 * beta;end 
	if lsiter<1, t0 = t0 / beta; end

    % update x
	x = (x + t*dx);
	% print some numbeZrs	
    if param.display,
        fprintf(' ite = %d \t cost = %f \n',k,f1);
    end

    %conjugate gradient calculation
	g1 = grad(x,param);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	% stopping criteria (to be improved)
	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end

end
return;

function res = objective(x,dx,t,param) %**********************************

% L2-norm part
w=param.E*(x+t*dx)-param.y;
L2Obj=w(:)'*w(:);

% L1-norm part FFT
if param.L1Weight
   w = param.W*(x+t*dx); 
   L1Obj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
else
    L1Obj=0;
end

% L1-norm part PCA
if param.L1Weight_pca 
   [w,~,~] = tempPCA(x+t*dx);
   L1Obj_pca = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
else
    L1Obj_pca=0;
end

% TV part
if param.TVWeight
   w = param.TV*(x+t*dx); 
   TVObj = sum((w(:).*conj(w(:))+param.l1Smooth).^(1/2));
else
    TVObj=0;
end
% objective function
res=L2Obj+param.L1Weight*L1Obj+param.TVWeight*TVObj+param.L1Weight_pca*L1Obj_pca;

function g = grad(x,param)%***********************************************

% L2-norm part
L2Grad = 2.*(param.E'*(param.E*x-param.y));

% L1-norm part FFT
if param.L1Weight
    w = param.W*x;
    L1Grad = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    L1Grad=0;
end


% L1-norm part PCA
if param.L1Weight_pca
    w = param.W*x;
    L1Grad_pca = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    L1Grad_pca=0;
end



% TV part
if param.TVWeight
    w = param.TV*x;
    TVGrad = param.TV'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
else
    TVGrad=0;
end

% complete gradient
g=L2Grad+param.L1Weight*L1Grad+param.TVWeight*TVGrad+param.L1Weight_pca*L1Grad_pca;

