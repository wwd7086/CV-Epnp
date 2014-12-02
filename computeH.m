function H2to1 = computeH(p1, p2) % p1, p2 should be n*2 matrix

n = size(p1,1);
if (n<4)
    H2to1 = 0;
    return
else
    T1 = zeros(2*n,n); % form a matrix to help generate A
    k = 0:1:n-1;
    T1(2*n*k+2*k+1) = 1;
    T1(2*n*k+2*k+2) = 1;
    A_right_2 = T1*p2;
    A_right_2 = cat(2, A_right_2, ones(2*n, 1));

    T2 = zeros(2*n,2*n);
    T2(4*n*k+2*k+1) = p1(k+1);
    T2(4*n*k+2*k+2*n+2) = p1(n+k+1);
    A_right = -T2*A_right_2;

    A_left = zeros(6,2*n);
    A_left(12*k+1) = p2(k+1);
    A_left(12*k+2) = p2(n+k+1);
    A_left(12*k+3) = 1;
    A_left(12*k+9+1) = p2(k+1);
    A_left(12*k+9+2) = p2(n+k+1);
    A_left(12*k+9+3) = 1;

    A=cat(2, A_left', A_right);
    

    A_0=A'*A;
    [V,D] = eig(A_0);
%         k2 = 0:1:2*n-1;
%         d(k2+1) = D(2*n*k2+k2+1);
%         [C, I] = min(d);
    target_v = V(:,1);
    H2to1 = reshape(target_v, 3, 3);
    H2to1 = H2to1';
    scale_H = abs(det(H2to1))^(1/3);
    if(det(H2to1)>0)
        H2to1 = bsxfun(@times, H2to1, 1/scale_H);
    else
        H2to1 = bsxfun(@times, H2to1, -1/scale_H);
    end
end


