function [inliers,M] = distpnp(M,X,t)

K = [961.9668 0 0; 0 966.7282 0; 660.5078 356.6476 1.0000]';
Xin = X(1:3,:);
Xout = X(4:end,:);

Xp = K * M * [Xin;ones(1,size(Xin,2))];
Xp = bsxfun(@rdivide, Xp, Xp(3,:));
Xp = Xp(1:2,:);

dist = sqrt(sum((Xout - Xp).^2));
inliers = find(dist<t);

end