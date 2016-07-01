function [ obj ] = getBNSst( dataA, l, delta, gA, numpool, opts)
%%
if nargin < 6
    opts = [];
    opts.niter = 20000;
    opts.br = 10000;
    opts.nu = 0;
    opts.lambda = 0;
    opts.eta1 = -0.5;
    opts.etaS = 0;
    opts.etaT = 0;
    opts.fixetaS = 0;
    opts.fixetaT = 0;
    opts.parallel_t=0;
    opts.parallel_s=0;
    opts.parallel_line=0;
end
%%
[nt, nbr] = size(dataA);
[~, d] = size(dataA{1,1});
bseA = ones(d, d, nt, nbr);
%%
for t = 1:nt 
    for br = 1:nbr
        data = dataA{t, br};
        [nrepli, ~] = size(data);
        if nrepli<=2     
        else
            for y = 1:d
                Y = data(:,y);
                bseA(y, [1:(y-1) (y+1):d], t ,br) = std(Y)*l;
            end
        end  
    end
end
%%
dd = zeros(d, d, nt, nbr);
nrepliA = zeros(nt, nbr);
for t = 1:nt 
    for br = 1:nbr
        data = dataA{t, br};
        [nrepliA(t, br), ~] = size(data);
        if nrepliA(t, br) > 2
            dd(:,:,t, br) = data'*data;
        end
    end
end
%%
if gA(1) == -1
    gA = binornd(1, 0.2, d, d, nt, nbr);
    %make it symmetric
    for t = 1:nt 
        for br = 1:nbr
            if nrepliA(t, br) > 2
                gtmp = gA(:, :, t, br); 
                gtmp = gtmp - tril(gtmp);
                gtmp = gtmp + gtmp';
                gA(:, :, t, br) = gtmp; 
            else
                gA(:, :, t, br) = 0.5; 
            end
        end
    end
end

%%
nu = opts.nu; lambda = opts.lambda;
bAnull = zeros(d, d, nt, nbr);
gAsum = bAnull;
count = 0;
eta1 = opts.eta1; etaS = opts.etaS; etaT = opts.etaT; 
etaSA = []; etaTA = [];
sA = ones(d, nt, nbr);
shapeA = repmat( (nu+nrepliA)/2, 1, 1, d);
shapeA = permute(shapeA, [3, 1, 2]);
scaleA = sA;
if opts.parallel_t == 1 || opts.parallel_s == 1 || opts.parallel_line==1
    parpool(numpool);
end
for iter = 1:opts.niter
    if mod(iter, 200)==0
        fprintf('  completed %d%% \r',round(iter/opts.niter*100));
    end
    %%update beta
    if opts.parallel_t==1
        parfor t = 1:nt 
            dd1 = squeeze(dd(:,:,t,:));
            for br = 1:nbr
                data = dataA{t, br};
                if nrepliA(t, br) > 2
                    for y = 1:d
                        Y = data(:,y);
                        g = gA(y, :, t, br);
                        g = g([1:(y-1) (y+1):d]);
                        bse = bseA(y,:,t, br);
                        bse = bse([1:(y-1) (y+1):d]);
                        s = sA(y, t, br);
                        lab = (delta^2).^abs(g-1);
                        dr = bse.^2 .* lab;
                        A = dd1([1:(y-1) (y+1):d],[1:(y-1) (y+1):d],br); 
                        A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                        R = chol(A);
                        z = (R')\(X'*Y);
                        bnulltmp = R\(s.*randn(d-1, 1)+z);
                        bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
                        bAnull(y, :, t, br)= bnulltmp;
                        scaleA(y, t, br) = 1/(nu*lambda/2 + sum(( Y - data*bnulltmp ).^2)/2);
                    end
                end
            end
        end
    elseif opts.parallel_s==1
        parfor br = 1:nbr 
            for t = 1:nt
                data = dataA{t, br};
                dd1 = squeeze(dd(:,:,:,br));
                if nrepliA(t, br) > 2
                    for y = 1:d
                        Y = data(:,y);
                        g = gA(y, :, t, br);
                        g = g([1:(y-1) (y+1):d]);
                        bse = bseA(y,:,t, br);
                        bse = bse([1:(y-1) (y+1):d]);
                        s = sA(y, t, br);
                        lab = (delta^2).^abs(g-1);
                        dr = bse.^2 .* lab;
                        A = dd1([1:(y-1) (y+1):d],[1:(y-1) (y+1):d],t); 
                        A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                        R = chol(A);
                        z = (R')\(X'*Y);
                        bnulltmp = R\(s.*randn(d-1, 1)+z);
                        bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
                        bAnull(y, :, t, br)= bnulltmp;
                        scaleA(y, t, br) = 1/(nu*lambda/2 + sum(( Y - data*bnulltmp ).^2)/2);
                    end
                end
            end
        end
    elseif opts.parallel_line==1
        for t = 1:nt 
            for br = 1:nbr
                data = dataA{t, br};
                dd1 = squeeze(dd(:,:,t,br));
                if nrepliA(t, br) > 2
                    parfor y = 1:d
                        Y = data(:,y);
                        g = gA(y, :, t, br);
                        g = g([1:(y-1) (y+1):d]);
                        bse = bseA(y,:,t, br);
                        bse = bse([1:(y-1) (y+1):d]);
                        s = sA(y, t, br);
                        lab = (delta^2).^abs(g-1);
                        dr = bse.^2 .* lab;            
                        A = dd1([1:(y-1) (y+1):d],[1:(y-1) (y+1):d]); 
                        A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                        R = chol(A);
                        z = (R')\(X'*Y);
                        bnulltmp = R\(s.*randn(d-1, 1)+z);
                        bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
                        bAnull(y, :, t, br)= bnulltmp;
                        scaleA(y, t, br) = 1/(nu*lambda/2 + sum(( Y - data*bnulltmp ).^2)/2);
                    end
                end
            end
        end
    else
        for t = 1:nt 
            for br = 1:nbr
                data = dataA{t, br};
                if nrepliA(t, br) > 2
                    for y = 1:d
                        Y = data(:,y);
                        X = data(:,[1:(y-1) (y+1):d]); 
                        g = gA(y, [1:(y-1) (y+1):d], t, br);
                        bse = bseA(y,[1:(y-1) (y+1):d], t, br);
                        s = sA(y, t, br);
                        lab = (delta^2).^abs(g-1);
                        dr = bse.^2 .* lab;
                        A = dd([1:(y-1) (y+1):d],[1:(y-1) (y+1):d],t,br); 
                        A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                        R = chol(A);
                        z = (R')\(X'*Y);
                        bAnulltmp = R\(s.*randn(d-1, 1)+z);
                        bAnull(y, [1:(y-1) (y+1):d], t, br)= bAnulltmp;
                        scaleA(y, t, br) = 1/(nu*lambda/2 + sum(( Y - X*bAnulltmp ).^2)/2);
                    end
                end
            end
        end
    end   
    
    %%update sigma
    sA = 1./sqrt(gamrnd(shapeA, scaleA));
    
    %%update latent state
    gA = updategASTmiss(gA, nrepliA, nt, nbr, d, bAnull, bseA, delta, eta1, etaS, etaT, 1);
    
    tmp = getxySTmiss(gA, nt, nbr, d, 1);
    %update etaS
    if opts.fixetaS~=1
      etaS = getetaSmhst(tmp.xs, tmp.xt, tmp.y, eta1, etaS, etaT); 
    end
    %update etaT
    if opts.fixetaT~=1
      etaT = getetaTmhst(tmp.xs, tmp.xt, tmp.y, eta1, etaS, etaT);
    end
    if iter >= opts.br 
      count = count + 1;
      gAsum = gAsum + gA;
      etaSA = [etaSA etaS];
      etaTA = [etaTA etaT];
    end
end
%%change the diagonal elements
for t = 1:nt
    for br = 1:nbr
        gtmp = gAsum(:, :, t, br);
        gtmp = gtmp - diag(diag(gtmp)) + eye(d)*count;
        gAsum(:, :, t, br) = gtmp;
    end
end


%%output
obj.tau1=bseA; 
obj.tau0=bseA*delta;
obj.postprob=gAsum/count; 
obj.etaSA=etaSA; 
obj.etaTA=etaTA;

end