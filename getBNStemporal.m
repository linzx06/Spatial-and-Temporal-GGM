function [ obj ] = getBNStemporal( dataA, l, delta, gA, opts)
%%
if nargin < 5
    opts = [];
    opts.niter = 20000;
    opts.br = 10000;
    opts.nu = 0;
    opts.lambda = 0;
    opts.eta1 = -0.5;
    opts.etaT = 0;
    opts.fixetaT = 0;
    opts.constrain = 1;
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
    if opts.constrain
        gA = binornd(1, 0.5, d, d, nsam);
        %make it symmetric
        for sam = 1:nsam 
            gtmp = gA(:, :, sam); 
            gtmp = gtmp - tril(gtmp);
            gtmp = gtmp + gtmp';
            gA(:, :, sam) = gtmp; 
        end
    else
        gA = binornd(1, 0.5, d, d, nsam);
    end
end
%%
AXX = zeros(d-1, d-1, d, nsam);
AXY = zeros(d-1, d, nsam);
nrepliA = zeros(nsam, 1);
for sam = 1:nsam 
    data = dataA{sam};
    [nrepliA(sam), ~] = size(data);
    for y = 1:d
        Y = data(:,y);
        X = data(:,[1:(y-1) (y+1):d]);
        AXX(:,:,y, sam) = X'*X;
        AXY(:,y, sam) = X'*Y;
    end
end
%%
nu = opts.nu; lambda = opts.lambda;
bAnull = zeros(d, d, nsam);
gAsum = bAnull;
count = 0;
eta1 = opts.eta1; etaT = opts.etaT; 
etaTA = [];
sA = ones(d, nsam);
shapeA = repmat((nu+nrepliA)/2, 1, d);
shapeA = permute(shapeA, [2, 1]);
scaleA = sA;
for iter = 1:opts.niter
    if mod(iter, 200)==0
        fprintf('  completed %d%% \r',round(iter/opts.niter*100));
    end
    %update beta
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
                A = AXX(:,:,y, sam); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\AXY(:, y, sam);
                bAnulltmp = R\(s.*randn(d-1, 1)+z);
                bAnull(y, [1:(y-1) (y+1):d], sam)= bAnulltmp;
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - X*bAnulltmp ).^2)/2);
            end
        end
    elseif opts.parallel_graph == 1 && opts.parallel_line == 0
        parfor sam = 1:nsam 
            data = dataA{sam};
            gAtmp = gA(:, :, sam);
            bseAtmp = bseA(:,:, sam);
            bAnulltmp = zeros(d, d);
            for y = 1:d
                Y = data(:,y);
                g = gAtmp(y, [1:(y-1) (y+1):d]);
                bse = bseAtmp(y,[1:(y-1) (y+1):d]);
                s = sA(y, sam);
                lab = (delta^2).^abs(g-1);
                dr = bse.^2 .* lab;
                A = AXX(:,:,y, sam); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\AXY(:, y, sam);
                bAnulltmp(y, [1:(y-1) (y+1):d])= R\(s.*randn(d-1, 1)+z);
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - data*(bAnulltmp(y, :))' ).^2)/2);
            end
            bAnull(:, :, sam) = bAnulltmp;
        end
    elseif opts.parallel_graph == 0 && opts.parallel_line == 1
        for sam = 1:nsam 
            data = dataA{sam};
            gAtmp = gA(:, :, sam);
            bseAtmp = bseA(:,:, sam);
            bAnulltmp = zeros(d, d);
            parfor y = 1:d
                Y = data(:,y);
                g = gAtmp(y, :);
                g = g([1:(y-1) (y+1):d]);
                bse = bseAtmp(y,:);
                bse = bse([1:(y-1) (y+1):d]);
                s = sA(y, sam);
                lab = (delta^2).^abs(g-1);
                dr = bse.^2 .* lab;
                A = AXX(:,:,y, sam); 
                A((1:d:((d-1)^2))) = A((1:d:((d-1)^2))) + s^2./dr;
                R = chol(A);
                z = (R')\AXY(:, y, sam);
                bnulltmp = R\(s.*randn(d-1, 1)+z);
                bnulltmp = [bnulltmp(1:(y-1)); 0 ;bnulltmp(y:(d-1))];
                bAnulltmp(y, :) = bnulltmp;
                scaleA(y, sam) = 1/(nu*lambda/2 + sum(( Y - data*(bAnulltmp(y, :))' ).^2)/2);
            end
            bAnull(:, :, sam) = bAnulltmp;
        end
    else
        error('parallel_graph and parallel_line cannot be both one');
    end
    
    %update sigma
    sA = 1./sqrt(gamrnd(shapeA, scaleA));
    
    %update the latent state
    gA = updategAT(gA, nsam, d, bAnull, bseA, delta, eta1, etaT, opts.constrain);
    
    tmp = getxyT(gA, nsam, d, opts.constrain);
    %update etaT
    if opts.fixetaT~=1
      etaT = getetaTmh(tmp.xt, tmp.y, eta1, etaT);
    end

    if iter >= opts.br && mod(iter, 5)==0
      count = count + 1;
      gAsum = gAsum + gA;
      etaTA = [etaTA etaT];
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
obj.etaTA=etaTA; 
end