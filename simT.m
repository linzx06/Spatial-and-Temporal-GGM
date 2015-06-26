function [ obj ] = simT(nrepli, ng, nt, prob, tran, value, opts)
%%
if nargin < 9
    opts = [];
    opts.minv = 0.1;
    opts.maxv = 0.4;
    opts.ratio = 1;
    opts.diagA = 0.5;
end
%%
%upper-diagonal elements, otherwise -1
supt = zeros(ng, ng, nt);
for t = 1:nt
    if t == 1
        tmp = binornd(1, prob, ng, ng);
        tmp = tmp - tril(tmp) - tril(ones(ng));
        supt(:, :, t) = tmp;
    else
        inx1 = find(supt(:, :, t-1)==1);
        inx0 = find(supt(:, :, t-1)==0);
        nc = round(length(inx1)*tran);
        
        inx1 = randsample(inx1, length(inx1));
        inx0 = [inx0; inx1(1:nc)];
        inx1 = inx1((nc+1):length(inx1));
        
        inx0 = randsample(inx0, length(inx0));
        inx1 = [inx1; inx0(1:nc)];
        %inx0 = inx0((nc+1):length(inx0));
        
        tmp = zeros(ng);
        tmp(inx1) = 1;
        tmp = tmp - tril(ones(ng));
        supt(:, :, t) = tmp;
    end
end
%%
%symmetrize
for t = 1:nt
    tmp = supt(:, :, t); 
    tmp = tmp - tril(tmp);
    tmp = tmp + tmp'; tmp = tmp + eye(ng);
    supt(:, :, t) = tmp;
end

%%
omega = zeros(ng, ng, nt);
if strcmp(value, 'different')
    for t = 1:nt
        vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
        vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
        tmp = vas .* supt(:, :, t); 
        tmp = tmp - tril(tmp);
        tmp = tmp + tmp';
        omega(:, :, t) = tmp;
    end 
elseif strcmp(value, 'same')
    vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
    vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
    for t = 1:nt
        tmp = vas .* supt(:, :, t); 
        tmp = tmp - tril(tmp);
        tmp = tmp + tmp';
        omega(:, :, t) = tmp;
    end
else
   error('Please input a valid value option, different or same') 
end
%%
%make it positive definite
for t = 1:nt
    omega(:, :, t) = omega(:, :, t) + diag(sum(abs(omega(:, :, t))) + opts.diagA);
end
%%
%check positve definiteness
for t = 1:nt
    if all(eig(omega(:, :, t)) > 0)
    else
        error('The matrix is not positive definite')
    end
end
%%
%simulate data
dataA = cell(nt, 1);
for t = 1:nt
	dataA{t} = mvnrnd(zeros(1, ng), inv(omega(:, :, t)), nrepli); 
end
%%output
obj.supA = supt;
obj.omega = omega;
obj.dataA = dataA;

end

