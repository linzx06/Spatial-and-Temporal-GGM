function [ obj ] = updategAS( gA, nsam, d, bAnull, bseA, delta, eta1, etaS, constrain )
%%
gAsum = sum(2*gA-1, 3);
%%
for sam = 1:nsam
    gsumtmp = gAsum - (2*gA(:, :, sam)-1);
    pr0 = 1./(exp(eta1 + etaS*gsumtmp) + 1);
    pr1 = 1 - pr0;
    if constrain~=1
        p1 = normpdf(bAnull(:, :, sam), 0, bseA(:, :, sam));
        p0 = normpdf(bAnull(:, :, sam), 0, bseA(:, :, sam)/delta);
        gA(:, :, sam) = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
    else 
        p1 = normpdf(bAnull(:, :, sam), 0, bseA(:, :, sam)) .* normpdf((bAnull(:, :, sam))', 0, (bseA(:, :, sam))');
        p0 = normpdf(bAnull(:, :, sam), 0, bseA(:, :, sam)*delta) .* normpdf((bAnull(:, :, sam))', 0, (bseA(:, :, sam))'*delta);
        gtmp = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
        gtmp = gtmp - tril(gtmp);
        gtmp = gtmp + gtmp';
        gA(:, :, sam) = gtmp;
    end
    gAsum = gsumtmp + (2*gA(:, :, sam)-1);
end
 
%%
obj = gA;

end
