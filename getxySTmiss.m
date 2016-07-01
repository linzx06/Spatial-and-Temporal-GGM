function [ obj ] = getxySTmiss( gA, nt, nbr, d, constrain )
%%
if constrain == 1
    gv = zeros(d*(d-1)/2, nt, nbr); xt = gv;
    for t = 1:nt
        for br = 1:nbr
            gtmp = gA(:, :, t, br);
            gtmp = gtmp - tril(gtmp) - tril(ones(d));
            gv(:, t, br) = gtmp(gtmp~=-1);
        end
    end
else
    gv = ones(d*(d-1), nt, nbr); xt = gv;
    for t = 1:nt
        for br = 1:nbr
            gtmp = gA(:, :, t, br);
            gtmp = gtmp - diag(diag(gtmp)) - eye(d);
            gv(:, t, br) = gtmp(gtmp~=-1);
        end
    end
end
%%
xs = repmat(sum(2*gv-1, 3), 1, 1, nbr) - (2*gv-1);
xt(:, 2:nt, :) = xt(:, 2:nt, :) + 2*gv(:, 1:(nt-1), :) - 1;
xt(:, 1:(nt-1), :) = xt(:, 1:(nt-1), :) + 2*gv(:, 2:nt, :) - 1;
%%
xs = reshape(xs, numel(xs), 1);
xt = reshape(xt, numel(xt), 1);
gv = reshape(gv, numel(gv), 1);
%%exclude the missing values
xs = xs(gv~=0.5);
xt = xt(gv~=0.5);
gv = gv(gv~=0.5);
%%
obj.xs = xs;
obj.xt = xt;
obj.y = gv;
end