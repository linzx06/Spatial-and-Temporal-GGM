function [ obj ] = simST(nrepli, ng, nt, nbr, prob, tran, purt, value, opts)
%%
if nargin < 9
    opts = [];
    opts.minv = 0.1;
    opts.maxv = 0.4;
    opts.add = 0.5;
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
%add perturbation and symmetrize
supA = zeros(ng, ng, nt, nbr);
for t = 1:nt
    tmpt = supt(:, :, t);  
    inx1 = find(tmpt==1);
    inx0 = find(tmpt==0);
    nc = round(length(inx1)*purt);
    for br = 1:nbr
        inx1tmp = randsample(inx1, length(inx1));
        inx0tmp = [inx0; inx1tmp(1:nc)];
        inx1tmp = inx1tmp((nc+1):length(inx1tmp));
        
        inx0tmp = randsample(inx0tmp, length(inx0tmp));
        inx1tmp = [inx1tmp; inx0tmp(1:nc)];
        tmp = zeros(ng); tmp(inx1tmp) = 1; 
        tmp = tmp + tmp'; tmp = tmp + eye(ng);
        supA(:, :, t, br) = tmp;
    end
end
%%
omega = zeros(ng, ng, nt, nbr);
if strcmp(value, 'different')
    for t = 1:nt
        for br = 1:nbr
           vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
           vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
           tmp = vas .* supA(:, :, t, br); 
           tmp = tmp - tril(tmp);
           tmp = tmp + tmp';
           omega(:, :, t, br) = tmp;
        end
    end 
elseif strcmp(value, 'same')
    vas = rand(ng)*(opts.maxv-opts.minv) + opts.minv;
    vas = vas .* (binornd(1,0.5,ng,ng)-0.5)*2; 
    for t = 1:nt
        for br = 1:nbr
            tmp = vas .* supA(:, :, t, br); 
            tmp = tmp - tril(tmp);
            tmp = tmp + tmp';
            omega(:, :, t, br) = tmp;
        end
    end
else
   error('Please input a valid value option, different or same') 
end
%%
%make it positive definite
for t = 1:nt
    for br = 1:nbr
        omega(:, :, t, br) = omega(:, :, t, br) + diag(sum(abs(omega(:, :, t, br))) + opts.add);
    end
end
%%
%check positve definiteness
for t = 1:nt
    for br = 1:nbr
        if all(eig(omega(:, :, t, br)) > 0)
        else
            error('The matrix is not positive definite')
        end
    end
end
%%
%simulate data
dataA = cell(nt, nbr);
for t = 1:nt
    for br = 1:nbr 
       dataA{t, br} = mvnrnd(zeros(1, ng), inv(omega(:, :, t, br)), nrepli); 
    end
end
%%output
obj.supA = supA;
obj.omega = omega;
obj.dataA = dataA;

end
