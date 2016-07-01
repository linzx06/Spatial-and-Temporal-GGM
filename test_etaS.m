cd /Users/zhixianglin/Documents/Graphical Models/revision_code/comparison_wishart

filenameD = 'd20_n100_change20.mat';
dataA = load(filenameD);
dataA = dataA.dataA;


%%%%BNS
l = 0.1; delta = 0.1;
opts = [];
opts.niter = 20000;
opts.br = 10000;
opts.nu = 0;
opts.lambda = 0;
opts.eta1 = -2;
opts.etaS = 0;
opts.fixetaS = 0;
opts.parallel_graph = 0;
opts.parallel_line = 0;

result = getBNSspatialnew( dataA.dataA, l, delta, -1, 4, opts);
gest = result.postprob;
d = 20;
n = 100;
nsam = 3;
prob = 0.1;
change = 0.2;
value = 'same';
filenameFP = strcat(num2str(d), '_', num2str(n), '_nsam_', num2str(nsam), '_prob_', num2str(prob), '_change_', num2str(change),'_value_', num2str(value), '_l_', num2str(l), '_delta_', num2str(delta), 'etaS0t1_FPb.csv');
filenameTP = strcat(num2str(d), '_', num2str(n), '_nsam_', num2str(nsam), '_prob_', num2str(prob), '_change_', num2str(change),'_value_', num2str(value), '_l_', num2str(l), '_delta_', num2str(delta), 'etaS0t1_TPb.csv');
supTrue = dataA.supA;
TPb = [];
FPb =[];
for cf = [1/1000: 1/1000:1]
	ztmp = (gest>=cf);
	TP = sum(supTrue(ztmp == 1)) - d*nsam;
	FP = sum((supTrue(ztmp == 1))==0);
	TPb = [TPb, TP];
	FPb = [FPb, FP];
end
dlmwrite(filenameFP, FPb, '-append');
dlmwrite(filenameTP, TPb, '-append');


filenameD = 'd20_n100_change100.mat';
dataA = load(filenameD);
dataA = dataA.dataA;

result1 = getBNSspatialnew( dataA.dataA, l, delta, -1, 4, opts);
gest = result1.postprob;
d = 20;
n = 100;
nsam = 3;
prob = 0.1;
change = 1;
value = 'same';
filenameFP = strcat(num2str(d), '_', num2str(n), '_nsam_', num2str(nsam), '_prob_', num2str(prob), '_change_', num2str(change),'_value_', num2str(value), '_l_', num2str(l), '_delta_', num2str(delta), 'etaS0t1_FPb.csv');
filenameTP = strcat(num2str(d), '_', num2str(n), '_nsam_', num2str(nsam), '_prob_', num2str(prob), '_change_', num2str(change),'_value_', num2str(value), '_l_', num2str(l), '_delta_', num2str(delta), 'etaS0t1_TPb.csv');
supTrue = dataA.supA;
TPb = [];
FPb =[];
for cf = [1/1000: 1/1000:1]
	ztmp = (gest>=cf);
	TP = sum(supTrue(ztmp == 1)) - d*nsam;
	FP = sum((supTrue(ztmp == 1))==0);
	TPb = [TPb, TP];
	FPb = [FPb, FP];
end
dlmwrite(filenameFP, FPb, '-append');
dlmwrite(filenameTP, TPb, '-append');