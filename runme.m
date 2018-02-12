clear all
load data
dat = ideal_filter(Craddock_results.timeseries',0.72,[0.01 0.1]);
%dat = ideal_filter(rand(size(Craddock_results.timeseries))',0.72,[0.01 0.1]);

h2 = angle(hilbert(dat(100:end-100,:)));
st = 6;
kk=0
for s = st:20
    s
    for k = 1:1000
        cluster = datasample(1:200,s,'Replace',false);
        phi = abs(metrics(h2(:,cluster)));
        meta(s-(st-1),k) = std(phi);
        sync(s-(st-1),k) = mean(phi);
        ent(s-(st-1),k) = entropy(phi);
        sz(s-(st-1),k) = s;
        ids{s,k} = cluster;
        kk = kk+1;
    end
end
figure;
hold on;
scatter(sync(:),ent(:),20,sz(:))
text(sync(:),ent(:),cellstr(num2str((1:kk)')));
