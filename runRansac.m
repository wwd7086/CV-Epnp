clear; close;
load('pointSet1');
load('cameraParams');

%5*d, first 3 rows is the 3D point in world frame
%last 2 rows is the 2D point in image
X = [fixedPoints';ones(1,size(fixedPoints,1));movingPoints'];

[M, inliers] = ransac(X,'fitpnp','distpnp','degenpnp',6,10,1);

M


