function [ obj ] = getBNSspatial( dataA, l, delta, gA, numpool, opts)
%%
if nargin < 6
    opts = [];
    opts.niter = 20000;
    opts.br = 10000;
    opts.nu = 0;
    opts.lambda = 0;
    opts.eta1 = -0.5;
    opts.etaS = 0;
    opts.fixetaS = 0;
    opts.parallel_graph = 0;
    opts.parallel_line = 0;
end
%%
nsam = length(dataA);
[~, d] = size(dataA{1});
%%
bseA = ones(d, d, nsam);
for i = 1:nsam
        data = dataA{i};
        for y = 1:d
            Y = data(:,y);      
            bseA(y, [1:(y-1) (y+1):d], i) = std(Y)*l;
        end
end
%%
if gA(1) == -1
    gA = binornd(1, 0.3, d, d, nsam);
    %make it symmetric
    for sam = 1:nsam 
        gtmp = gA(:, :, sam); 
        gtmp = gtmp - tril(gtmp);
        gtmp = gtmp + gtmp';
        gA(:, :, sam) = gtmp; 
    end
end
%%
dd = zeros(d, d, nsam);
nrepliA = zeros(nsam, 1);
for sam = 1:nsam 
    data = dataA{sam};
    [nrepliA(sam), ~] = size(data);
    dd(:,:,sam) = data'*data;
end
%%
nu = opts.nu; lambda = opts.lambda;
bAnull = zeros(d, d, nsam);
gAsum = bAnull;
count = 0;
eta1 = opts.eta1; etaS = opts.etaS; 
etaSA = [];
sA = ones(d, nsam);
shapeA = repmat((nu+nrepliA)/2, 1, d);
shapeA = permute(shapeA, [2, 1]);
scaleA = sA;
if opts.parallel_graph == 1 || opts.parallel_line == 1
    parpool(numpool);
end
for iter = 1:opts.niter
    if mod(iter, 200)==0
        fprintf('completed %d%% \r',round(iter/opts.niter*100));
    end
    %%update beta
    if opts.parallel_graph == 0 && opts.parallel_line == 0
        for sam = 1:nsam 
            data = dataA{sam};
            for y = 1:d
                Y = data(:,y);
                X = data(:,[1:(y-1) (y+1):d]); 
                g = gA(y, [1:(y-1) (y+1):d], sam);
                bse = bseA(y,[1:(y-1) (y+1):d], sam);
                s = sA(y, sam);
                lab = (delta^2).^abs(g-1);
                dr = bse.^2 .* lab;
                A = dd([1:(y-1) (y+1):d],[1:(y-1) (y+1):d],sam); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\(X'*Y);
                bAnulltmp = R\(s.*randn(d-1, 1)+z);
                bAnull(y, [1:(y-1) (y+1):d], sam)= bAnulltmp;
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - X*bAnulltmp ).^2)/2);
            end
        end
    elseif opts.parallel_graph == 1 && opts.parallel_line == 0
        parfor sam = 1:nsam 
            data = dataA{sam};
            dd1 = dd(:,:,sam);
            gAtmp = gA(:, :, sam);
            bseAtmp = bseA(:,:, sam);
            bAnulltmp = zeros(d, d);
            for y = 1:d
                X = data(:,[1:(y-1) (y+1):d]); 
                Y = data(:,y);
                g = gAtmp(y, [1:(y-1) (y+1):d]);
                bse = bseAtmp(y,[1:(y-1) (y+1):d]);
                s = sA(y, sam);
                lab = (delta^2).^abs(g-1);
                dr = bse.^2 .* lab;
                A = dd1([1:(y-1) (y+1):d],[1:(y-1) (y+1):d]); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\(X'*Y);
                bAnulltmp(y, [1:(y-1) (y+1):d])= R\(s.*randn(d-1, 1)+z);
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - data*(bAnulltmp(y, :))' ).^2)/2);
            end
            bAnull(:, :, sam) = bAnulltmp;
        end
    elseif opts.parallel_graph == 0 && opts.parallel_line == 1
        for sam = 1:nsam 
            data = dataA{sam};
            dd1 = dd(:,:,sam);
            gAtmp = gA(:, :, sam);
            bseAtmp = bseA(:,:, sam);
            bAnulltmp = zeros(d, d);
            parfor y = 1:d
                X = data(:,[1:(y-1) (y+1):d]); 
                Y = data(:,y);
                g = gAtmp(y, :);
                g = g([1:(y-1) (y+1):d]);
                bse = bseAtmp(y,:);
                bse = bse([1:(y-1) (y+1):d]);
                s = sA(y, sam);
                lab = (delta^2).^abs(g-1);
                dr = bse.^2 .* lab;
                A = dd1([1:(y-1) (y+1):d],[1:(y-1) (y+1):d]); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\(X'*Y);
                bnulltmp = R\(s.*randn(d-1, 1)+z);
                bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
                bAnulltmp(y, :) = bnulltmp;
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - data*bnulltmp ).^2)/2);
            end
            bAnull(:, :, sam) = bAnulltmp;
        end
    else
        error('parallel_graph and parallel_line cannot be both one');
    end
    
    %update sigma
    sA = 1./sqrt(gamrnd(shapeA, scaleA));
    
    %update the latent state
    gA = updategAS(gA, nsam, d, bAnull, bseA, delta, eta1, etaS, 1);
    
    tmp = getxyS(gA, nsam, d, 1);
    %update etaS
    if opts.fixetaS~=1
      etaS = getetaSmh(tmp.xs, tmp.y, eta1, etaS);
    end
    if iter >= opts.br 
      count = count + 1;
      gAsum = gAsum + gA;
      etaSA = [etaSA etaS];  
    end
end
%%change the diagonal elements
for sam = 1:nsam
    gtmp = gAsum(:, :, sam);
    gtmp = gtmp - diag(diag(gtmp)) + eye(d)*count;
    gAsum(:, :, sam) = gtmp;
end

%%output
obj.tau1=bseA; 
obj.tau0=bseA*delta;
obj.postprob=gAsum/count;  
obj.etaSA=etaSA; 
end