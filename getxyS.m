function [ obj ] = getxyS( gA, nsam, d, constrain )
%%
if constrain == 1
    gv = zeros(d*(d-1)/2, nsam);
    for sam = 1:nsam
        gtmp = gA(:, :, sam);
        gtmp = gtmp - tril(gtmp) - tril(ones(d));
        gv(:, sam) = gtmp(gtmp~=-1);
    end
else
    gv = ones(d*(d-1), nsam); 
    for sam = 1:nsam
        gtmp = gA(:, :, sam);
        gtmp = gtmp - diag(diag(gtmp)) - eye(d);
        gv(:, sam) = gtmp(gtmp~=-1);
    end
end
%%
xs = repmat(sum(2*gv-1, 2), 1, nsam) - (2*gv-1);
%%
xs = reshape(xs, numel(xs), 1);
gv = reshape(gv, numel(gv), 1);
%%
obj.xs = xs;
obj.y = gv;
end

