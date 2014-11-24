clear; close;
load('pointSet1');

%Camera Intrinstic parameter
K = [1 0 0; 0 1 0; 0 0 1];

%All points in image
%Pi = [ 1 4 7 6 3 2; 3 4 8 2 5 6];
Pi = movingPoints';

%All points in world frame
%Pw = [1 3 5 9 6 7; 5 8 6 4 3 2; 6 4 1 3 2 4]; %3*6
Pw = [fixedPoints';ones(1,size(fixedPoints,1))];

[R,T] = Epnp(K,Pi,Pw);

tPw = [Pw;ones(1,size(Pw,2))];
tPi = K*[R T]*tPw;

Pi
tPi