function [ obj ] = simsig(nrepli, d, prob, opts)
%%the 
if nargin < 4
    opts = [];
    opts.maxv = 0.4;
    opts.minv = 0.1;
    opts.add = 0.5;
end
%%
vas = rand(d)*(opts.maxv-opts.minv) + opts.minv;
vas = vas .* (binornd(1,0.5,d,d)-0.5)*2;
sup = binornd(1,prob,d,d);
omega = vas .* sup;
omega = omega - tril(omega,-1);
omega = omega + omega';
[nRows,nCols] = size(omega);
omega(1:(nRows+1):nRows*nCols) = 1;

%%
posdef = 0;
while posdef~=1
    if (sum(eig(omega)<=0));
        omega = omega + eye(d).* 0.01;
    else 
        posdef=1;
    end
end
omega = omega + eye(d) .* opts.add;

%%
obj.data = mvnrnd(zeros(1, d), inv(omega), nrepli);
obj.omega = omega;
obj.sup = (omega~=0);

end
