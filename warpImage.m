function [H2to1,warpedImg] = warpImage(img1, img2, pts)

p1 = pts(:,1:2);
p2 = pts(:,3:4);
H2to1 = computeH(p1, p2);
save('H.mat', 'H2to1');

p2_homo = cat(2, p2, ones(size(p2,1),1));
p2_trans_homo = H2to1*p2_homo';
p2_trans = p2_trans_homo(1:2,:);
p2_trans(1,:) = bsxfun(@rdivide, p2_trans(1,:), p2_trans_homo(3,:));
p2_trans(2,:) = bsxfun(@rdivide, p2_trans(2,:), p2_trans_homo(3,:));

err = bsxfun(@minus, p1', p2_trans);
rmse = sqrt(sum(sum(bsxfun(@power, err, 2)))/size(err,2)); % root mean squared err

outSize = [size(img1,1),size(img1,2)];
warpedImg = warpH(img2, H2to1, outSize, 0);
% imwrite(warpedImg, 'warped.jpg', 'jpg');

