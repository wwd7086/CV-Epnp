function M = fitpnp(X)
%M = [R T]
%X = 5*d  1:3 in   3:5 out
K = [961.9668 0 0; 0 966.7282 0; 660.5078 356.6476 1.0000]';
Xin = X(1:3,:);
Xout = X(4:end,:);

[R, T] = Epnp(K, Xout, Xin);

M = [R T];

end