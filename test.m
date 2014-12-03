clear;close;
addpath('OpenSURF_version1c');
addpath(genpath('MatchPoint'));
%-------- read reference image ------------------
I = imread('reference.jpg');
I2 = imread('view4.jpg');
%-------- make surface masks --------------------
[~,~, front_mask] = imread('front mask.png','PNG');
front_mask = bsxfun(@gt,front_mask,0);
[~,~, side_mask] = imread('side mask.png','PNG');
side_mask = bsxfun(@gt,side_mask,0);
[~,~, top_mask] = imread('top mask.png','PNG');
top_mask = bsxfun(@gt,top_mask,0);
masks{1}=front_mask;
masks{2}=side_mask;
masks{3}=top_mask;

%-------- load corner points --------------------
%load('cornerPointsSample');
%load('cornerPointsTest');
load('cornerPointsFinal');

%-------- calculate H matrixs for getting 3D coordinates ------------------
I_warped = zeros(200,200,3);
H2to1 = cell(1,3);
warpedImg = cell(1,3);
[H2to1{1},warpedImg{1}] = warpImage(I_warped, I, [image2d, front2d]); % image2d/front2d are all in [X,Y]
[H2to1{2},warpedImg{2}] = warpImage(I_warped, I, [image2d, side2d]);
[H2to1{3},warpedImg{3}] = warpImage(I_warped, I, [image2d, top2d]);

%-------- find the matching point --------------------
%[pts1, pts2] = findMatch(I,I2,0.0001,300);
%[pts1,pts2] = findMatchHaris(I,I2);
%save('matchingPoint1','pts1','pts2');

load('matchingclear');
subplot(1,5,1);
image(I);
hold on;

%------- plot the box -----------------------------------------
plot(front2d([1,2,3,4,1],1),front2d([1,2,3,4,1],2), 'c','LineWidth',3);
plot(side2d([1,2,3,4,1],1),side2d([1,2,3,4,1],2), 'c','LineWidth',3);
plot(top2d([1,2,3,4,1],1),top2d([1,2,3,4,1],2), 'c','LineWidth',3);

%-------- find the 3D coordinates of the feature point --------------------
objectSize = [200,200,200];
[x3d, surfaceNum] = find3dCoord(masks, H2to1, pts1, objectSize);
subplot(1,5,2);

%only keep the points on the surfaces
% pts1 -> x3d in the referene image 
% pts2 in the input image
x3d = x3d(:,surfaceNum~=0);
pts2 = pts2(:,surfaceNum~=0);
pts1 = pts1(:,surfaceNum~=0);
surfaceNum = surfaceNum(surfaceNum~=0);

%------- Visualize the 3D point --------------------------------
hold on;
image(warpedImg{1});
    axis([0 200 0 200]);
    axis ij;
subplot(1,5,3);
image(warpedImg{2});
    axis([0 200 0 200]);
    axis ij;
subplot(1,5,4);
image(warpedImg{3});
    axis([0 200 0 200]);
    axis ij;

plotpoints(1,pts1,'r',surfaceNum,x3d,'g');

%-------- calculate the R t -------------------------------
X = [x3d;pts2];
%[M, inliers] = ransac(X,'fitpnp','distpnp','degenpnp',6,100,1);
[R,t] = Epnp(getIntrinsic(),pts2,x3d);
M = [R t];
inliers = 1:13;

%plot the inliers on 2D image
pinlier = X(4:5,inliers);
p3inlier = X(1:3,inliers);
surfinlier = surfaceNum(inliers);

subplot(1,5,5);
image(I2);
hold on;
plotpoints(5,pinlier,'b',surfinlier,p3inlier,'b');

%------- Calculate reprojection -------------------------------
K = getIntrinsic();
R = M(1:3,1:3); % Rotation Matrix
t = M(1:3, 4); % Translation Vector

front3d = ones(4,4);
side3d = ones(4,4);
top3d = ones(4,4);

[front3d(1:3,:), ~] = find3dCoord(masks, H2to1, round(front2d'), objectSize);
[side3d(1:3,:), ~] = find3dCoord(masks, H2to1, round(side2d'), objectSize);
[top3d(1:3,:), ~] = find3dCoord(masks, H2to1, round(top2d'), objectSize);

front2dnew = K*[R t]*front3d; % in homography coordinates
side2dnew = K*[R t]*side3d;
top2dnew = K*[R t]*top3d;

front2dnew = bsxfun(@rdivide,front2dnew(1:2,:),front2dnew(3,:)); % in image coordinates
side2dnew = bsxfun(@rdivide,side2dnew(1:2,:),side2dnew(3,:));
top2dnew = bsxfun(@rdivide,top2dnew(1:2,:),top2dnew(3,:));

front2dnew = front2dnew';
side2dnew = side2dnew';
top2dnew = top2dnew';
%------- plot the new box -----------------------------------------
subplot(1,5,5);
hold on;

plot(front2dnew([1,2,3,4,1],1),front2dnew([1,2,3,4,1],2), 'c','LineWidth',3);
plot(side2dnew([1,2,3,4,1],1),side2dnew([1,2,3,4,1],2), 'c','LineWidth',3);
plot(top2dnew([1,2,3,4,1],1),top2dnew([1,2,3,4,1],2), 'c','LineWidth',3);

