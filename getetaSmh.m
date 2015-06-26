function [ obj ] = getetaSmh(xs, y, eta1, etaS, opts)
%%
if nargin < 5
    opts = [];
    opts.step = 0.05;
    opts.low = 0;
    opts.high = 2;
end
%%
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

lprobo = sum(y.*(eta1+etaS*xs) - log(1+exp(eta1+etaS*xs)));
lprobn = sum(y.*(eta1+etaSn*xs) - log(1+exp(eta1+etaSn*xs)));
prob = nto/otn*exp(lprobn-lprobo);
if prob>=randn(1)
    obj = etaSn;
else 
    obj = etaS;
end

end
