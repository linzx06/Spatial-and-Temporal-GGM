function [ obj ] = getxyT( gA, nt, d, constrain )
%%
if constrain == 1
    gv = zeros(d*(d-1)/2, nt); xt = gv;
    for t = 1:nt
        gtmp = gA(:, :, t);
        gtmp = gtmp - tril(gtmp) - tril(ones(d));
        gv(:, t) = gtmp(gtmp~=-1);
    end
else
    gv = ones(d*(d-1), nt, nbr); xt = gv;
    for t = 1:nt
        gtmp = gA(:, :, t);
        gtmp = gtmp - diag(diag(gtmp)) - eye(d);
        gv(:, t) = gtmp(gtmp~=-1);
    end
end
%%
xt(:, 2:nt) = xt(:, 2:nt) + 2*gv(:, 1:(nt-1)) - 1;
xt(:, 1:(nt-1)) = xt(:, 1:(nt-1)) + 2*gv(:, 2:nt) - 1;
%%
xt = reshape(xt, numel(xt), 1);
gv = reshape(gv, numel(gv), 1);
%%
obj.xt = xt;
obj.y = gv;

end

