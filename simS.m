function [ obj ] = simS(nrepli, ng, nsam, prob, change, value, opts)
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
supt = zeros(ng, ng, nsam);
for sam = 1:nsam
    if sam == 1
        tmp = binornd(1, prob, ng, ng);
        tmp = tmp - tril(tmp) - tril(ones(ng));
        supt(:, :, sam) = tmp;
    else
        inx1 = find(supt(:, :, 1)==1);
        inx0 = find(supt(:, :, 1)==0);
        nc = round(length(inx1)*change);
        
        inx1 = randsample(inx1, length(inx1));
        inx0 = [inx0; inx1(1:nc)];
        inx1 = inx1((nc+1):length(inx1));
        
        inx0 = randsample(inx0, length(inx0));
        inx1 = [inx1; inx0(1:nc)];
        %inx0 = inx0((nc+1):length(inx0));
        
        tmp = zeros(ng);
        tmp(inx1) = 1;
        tmp = tmp - tril(ones(ng));
        supt(:, :, sam) = tmp;
    end
end
%%
%symmetrize
for sam = 1:nsam
    tmp = supt(:, :, sam); 
    tmp = tmp - tril(tmp);
    tmp = tmp + tmp'; tmp = tmp + eye(ng);
    supt(:, :, sam) = tmp;
end

%%
omega = zeros(ng, ng, nsam);
if strcmp(value, 'different')
    for sam = 1:nsam
        vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
        vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
        tmp = vas .* supt(:, :, sam); 
        tmp = tmp - tril(tmp);
        tmp = tmp + tmp';
        omega(:, :, sam) = tmp;
    end 
elseif strcmp(value, 'same')
    vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
    vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
    for sam = 1:nsam
        tmp = vas .* supt(:, :, sam); 
        tmp = tmp - tril(tmp);
        tmp = tmp + tmp';
        omega(:, :, sam) = tmp;
    end
else
   error('Please input a valid value option, different or same') 
end
%%
%make it positive definite
for sam = 1:nsam
    omega(:, :, sam) = omega(:, :, sam) + diag(sum(abs(omega(:, :, sam))) + opts.diagA);
end
%%
%check positve definiteness
for sam = 1:nsam
    if all(eig(omega(:, :, sam)) > 0)
    else
        error('The matrix is not positive definite')
    end
end
%%
%simulate data
dataA = cell(nsam, 1);
for sam = 1:nsam
	dataA{sam} = mvnrnd(zeros(1, ng), inv(omega(:, :, sam)), nrepli); 
end
%%output
obj.supA = supt;
obj.omega = omega;
obj.dataA = dataA;

end

