clear; close;
load('matchingclear');
load('point3dclean');

addpath('OpenSURF_version1c');
addpath(genpath('MatchPoint'));
%Camera Intrinstic parameter
%K = cameraParams.IntrinsicMatrix';
%K = eye(3);
K = getIntrinsic();

%All points in image
%Pi = [ 1 4 7 6 3 2; 3 4 8 2 5 6];
Pi = pts2;
Pw = x3d;

[R,T] = Epnp(K,Pi,Pw);

tPi = K * [R T] * [Pw;ones(1,size(Pw,2))];
tPi = bsxfun(@rdivide, tPi, tPi(3,:));

%Pi
%tPi
R
T