function [ obj ] = updategASTmiss( gA, nrepliA, nt, nbr, d, bAnull, bseA, delta, eta1, etaS, etaT, constrain )
%%
gAsum = sum(2*gA-1,4);
%%
for t = 1:nt
    if t==1
        for br = 1:nbr
            if nrepliA(t, br) > 2
                gsumtmp = gAsum(:,:,t) - (2*gA(:, :, t, br)-1);
                pr0 = 1./(exp(eta1 + etaS*gsumtmp + etaT*(2*gA(:,:,t+1, br)-1)) + 1);
                pr1 = 1 - pr0;
                if constrain~=1
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br));
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta);
                    gA(:, :, t, br) = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                else 
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))');
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))'*delta);
                    gtmp = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                    gtmp = gtmp - tril(gtmp);
                    gtmp = gtmp + gtmp';
                    gA(:, :, t, br) = gtmp;
                end
                gAsum(:,:,t) = gsumtmp + (2*gA(:, :, t, br)-1);
            end
        end
    elseif t==nt
        for br = 1:nbr
            if nrepliA(t, br) > 2
                gsumtmp = gAsum(:,:,t) - (2*gA(:, :, t, br)-1);
                pr0 = 1./(exp(eta1 + etaS*gsumtmp + etaT*(2*gA(:,:,t-1, br)-1)) + 1);
                pr1 = 1 - pr0;
                if constrain~=1
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br));
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta);
                    gA(:, :, t, br) = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                else 
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))');
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))'*delta);
                    gtmp = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                    gtmp = gtmp - tril(gtmp);
                    gtmp = gtmp + gtmp';
                    gA(:, :, t, br) = gtmp;
                end
                gAsum(:,:,t) = gsumtmp + (2*gA(:, :, t, br)-1);
            end
        end
    else
        for br = 1:nbr
            if nrepliA(t, br) > 2
                gsumtmp = gAsum(:,:,t) - (2*gA(:, :, t, br)-1);
                pr0 = 1./(exp(eta1 + etaS*gsumtmp + etaT*(2*gA(:,:,t+1, br)-1+2*gA(:,:,t-1, br)-1)) + 1);
                pr1 = 1 - pr0;
                if constrain~=1
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br));
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta);
                    gA(:, :, t, br) = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                else 
                    p1 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))');
                    p0 = normpdf(bAnull(:,:,t, br), 0, bseA(:,:,t, br)*delta) .* normpdf((bAnull(:,:,t, br))', 0, (bseA(:,:,t, br))'*delta);
                    gtmp = (p1.*pr1./(p1.*pr1 + p0.*pr0)) >= rand(d);
                    gtmp = gtmp - tril(gtmp);
                    gtmp = gtmp + gtmp';
                    gA(:, :, t, br) = gtmp;
                end
                gAsum(:,:,t) = gsumtmp + (2*gA(:, :, t, br)-1);
            end
        end
    end   
end
%%
obj = gA;

end

