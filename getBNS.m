function [ obj ] = getBNS( data, q, l, delta, opts )
%%
if nargin < 5
    opts = [];
    opts.niter = 20000;
    opts.br = 10000;
    opts.nu = 0;
    opts.lambda = 0;
    opts.parallel = 0;
end

[nrepli, d] = size(data);
%%
bseA = ones(d);
for y = 1:d
    Y = data(:,y);
    bseA(y, [1:(y-1) (y+1):d]) = std(Y)*l;
end

%%
gA = ones(d);
gAsum = zeros(d);
sA = ones(1, d);
count = 0;
AXX = zeros(d-1, d-1, d);
AXY = zeros(d-1, d);
for y = 1:d
    Y = data(:,y);
    X = data(:,[1:(y-1) (y+1):d]);
    AXX(:,:,y) = X'*X;
    AXY(:,y) = X'*Y;
end;

if opts.parallel == 0
    for iter = 1:opts.niter
        if mod(iter, 200)==0
            fprintf('  completed %d%% \r',round(iter/opts.niter*100));
        end
        bnull = zeros(d);
        for y = 1:d
            Y = data(:,y);
            X = data(:,[1:(y-1) (y+1):d]); 
            g = gA(y, [1:(y-1) (y+1):d]);
            bse = bseA(y,[1:(y-1) (y+1):d]);
            s = sA(y);
            %update beta
            lab = (delta^2).^abs(g-1);
            dr = bse.^2 .* lab;
            A = AXX(:,:,y); 
            A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
            R = chol(A);
            z = (R')\AXY(:,y);
            bnulltmp = R\(s.*randn(d-1, 1)+z);
            bnull(y, [1:(y-1) (y+1):d]) = bnulltmp;
            %update sigma
            sA(y) = 1/sqrt(gamrnd((opts.nu+nrepli)/2, 1/(opts.nu*opts.lambda/2 + sum(( Y - X*bnulltmp ).^2)/2), 1));  
        end
        %update latent state
        p1 = normpdf(bnull, 0, bseA).*normpdf(bnull', 0, bseA')*q;
        p0 = normpdf(bnull, 0, bseA*delta).*normpdf(bnull', 0, bseA'*delta)*(1-q);
        gA = (p1./(p1 + p0)) >= rand(d);
        gA = gA - tril(gA);
        gA = gA + gA';
        if iter >= opts.br && mod(iter, 5)==0
            count = count + 1;
            gAsum = gAsum + gA;
        end
    end
else
   nu = opts.nu; lambda = opts.lambda;
   for iter = 1:opts.niter
        if mod(iter, 200)==0
            fprintf('  completed %d%% \r',round(iter/opts.niter*100));
        end
        bnull = zeros(d);
        parfor y = 1:d
            Y = data(:,y);
            g = gA(y, :);
            g = g([1:(y-1) (y+1):d]);
            bse = bseA(y,:);
            bse = bse([1:(y-1) (y+1):d]);
            s = sA(y);
            %update beta
            lab = (delta^2).^abs(g-1);
            dr = bse.^2 .* lab;
            A = AXX(:,:,y); 
            A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
            R = chol(A);
            z = (R')\AXY(:,y);
            bnulltmp = R\(s.*randn(d-1, 1)+z);
            bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
            bnull(:, y) = bnulltmp;
            %update sigma       
            sA(y) = 1/sqrt(gamrnd((nu+nrepli)/2, 1/(nu*lambda/2 + sum(( Y - data*bnulltmp ).^2)/2), 1));              
        end
        %update latent state
        p1 = normpdf(bnull, 0, bseA).*normpdf(bnull', 0, bseA')*q;
        p0 = normpdf(bnull, 0, bseA*delta).*normpdf(bnull', 0, bseA'*delta)*(1-q);
        gA = (p1./(p1 + p0)) >= rand(d);
        gA = gA - tril(gA);
        gA = gA + gA';
        if iter >= opts.br && mod(iter, 5)==0
            count = count + 1;
            gAsum = gAsum + gA;
        end
    end 
end
gAsum = gAsum - diag(diag(gAsum)) + eye(d)*count;

%%
obj.tau1=bseA; 
obj.tau0=bseA*delta;
obj.postprob=gAsum/count; 

end
