m = 20;
n = 10;
data = rand(m,n);

k = 3;
[idxs, mrlns] = multirelation(data, k, 'zvi');

