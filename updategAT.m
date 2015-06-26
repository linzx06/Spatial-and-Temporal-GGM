function [ obj ] = updategAT( gA, nt, d, bAnull, bseA, delta, eta1, etaT, constrain )
%%
for t = 1:nt
    if t == 1
        pr0 = 1./(exp(eta1 + etaT*(2*gA(:, :, t+1)-1)) + 1);
    elseif t == nt
        pr0 = 1./(exp(eta1 + etaT*(2*gA(:, :, t-1)-1)) + 1);
    else
        pr0 = 1./(exp(eta1 + etaT*(2*gA(:, :, t+1)-1+2*gA(:, :, t-1)-1)) + 1);
    end
    pr1 = 1 - pr0;
    if constrain~=1
        p1 = normpdf(bAnull(:, :, t), 0, bseA(:, :, t));
        p0 = normpdf(bAnull(:, :, t), 0, bseA(:, :, t)*delta);
        gA(:, :, t) = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
    else 
        p1 = normpdf(bAnull(:, :, t), 0, bseA(:, :, t)) .* normpdf((bAnull(:, :, t))', 0, (bseA(:, :, t))');
        p0 = normpdf(bAnull(:, :, t), 0, bseA(:, :, t)*delta) .* normpdf((bAnull(:, :, t))', 0, (bseA(:, :, t))'*delta);
        gtmp = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
        gtmp = gtmp - tril(gtmp);
        gtmp = gtmp + gtmp';
        gA(:, :, t) = gtmp;
    end
end
 
%%
obj = gA;

end
