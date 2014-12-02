function [inliers,M] = distpnp(M,X,t)

K = getIntrinsic();
Xin = X(1:3,:);
Xout = X(4:end,:);

Xp = K * M * [Xin;ones(1,size(Xin,2))];
Xp = bsxfun(@rdivide, Xp, Xp(3,:));
Xp = Xp(1:2,:);

dist = sqrt(sum((Xout - Xp).^2));
inliers = find(dist<t);

end