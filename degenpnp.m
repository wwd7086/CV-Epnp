function r = degenpnp(X)
% 0->is not degnerated
% 1->is degnerated
% X 5*N(sample)->6

Xin = X(1:3,:);
A = [Xin(:,1:3)',ones(3,1)];
L = null(A); %4*1

test = [Xin(:,4:end)',ones(size(X,2)-3,1)]; %N*4
err = mean(test*L);

if abs(err)>5
    r=0;
else
    r=1;
end

end