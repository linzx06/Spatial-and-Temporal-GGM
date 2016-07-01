function [ obj ] = getetaTmhst(xs, xt, y, eta1, etaS, etaT, opts)
%%
if nargin < 7
    opts = [];
    opts.step = 0.4;
    opts.low = 0;
    opts.high = 2;
end
%%
pd = makedist('Normal');
pd.mu = etaT;
pd.sigma = opts.step;
t = truncate(pd,opts.low,opts.high);
etaTn = random(t);
pd.mu = etaTn;
pd.sigma = opts.step;
tn = truncate(pd,opts.low,opts.high);
nto = pdf(tn, etaT);
otn = pdf(t, etaTn);

lprobo = sum(y.*(eta1+etaS*xs+etaT*xt) - log(1+exp(eta1+etaS*xs+etaT*xt)));
lprobn = sum(y.*(eta1+etaS*xs+etaTn*xt) - log(1+exp(eta1+etaS*xs+etaTn*xt)));
%prob = nto/otn*exp(lprobn-lprobo)*normpdf(btn, opts.m, opts.sd)/normpdf(bt, opts.m, opts.sd);
prob = nto/otn*exp(lprobn-lprobo);
if prob>=randn(1)
    obj = etaTn;
else 
    obj = etaT;
end

end


