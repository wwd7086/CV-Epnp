function M = fitpnp(X)
%M = [R T]
%X = 5*d  1:3 in   3:5 out
K = getIntrinsic();
Xin = X(1:3,:);
Xout = X(4:end,:);

[R, T] = Epnp(K, Xout, Xin);

M = [R T];

end