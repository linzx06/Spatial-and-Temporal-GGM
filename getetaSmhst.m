function [ obj ] = getetaSmhst(xs, xt, y, eta1, etaS, etaT, opts)
%%
if nargin < 7
    opts = [];
    opts.step = 0.05;
    opts.low = 0;
    opts.high = 2;
end

pd = makedist('Normal');
pd.mu = etaS;
pd.sigma = opts.step;
t = truncate(pd,opts.low,opts.high);
etaSn = random(t);
pd.mu = etaSn;
pd.sigma = opts.step;
tn = truncate(pd,opts.low,opts.high);
nto = pdf(tn, etaS);
otn = pdf(t, etaSn);

lprobo = sum(y.*(eta1+etaS*xs+etaT*xt) - log(1+exp(eta1+etaS*xs+etaT*xt)));
lprobn = sum(y.*(eta1+etaSn*xs+etaT*xt) - log(1+exp(eta1+etaSn*xs+etaT*xt)));
%prob = nto/otn*exp(lprobn-lprobo)*normpdf(bsn, opts.m, opts.sd)/normpdf(bs, opts.m, opts.sd);
prob = nto/otn*exp(lprobn-lprobo);
if prob>=randn(1)
    obj = etaSn;
else 
    obj = etaS;
end

end

