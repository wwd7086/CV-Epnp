function [R,T] = Epnp(K,Pi,Pw)

au = K(1);
av = K(5);
at = K(4);
u0 = K(7);
v0 = K(8);

numPw = size(Pw,2);

%pick 4 random points as control point
randCi = randperm(numPw);
randCi = randCi(1:4);
Cw = Pw(:,randCi);

%express all the point using control point
Acp = zeros(4,numPw);
for i=1:numPw

    acp = Cw\Pw(:,i);
    Acp(:,i) = acp;
    
end

%solve the control point in camer coordinate
M = zeros(numPw*2,12);
for i = 1:numPw
    acp = Acp(:,i);
    pi = Pi(:,i);
    
    m1c = [au, at, u0-pi(1)];
    m1 = [acp(1).*m1c, acp(2).*m1c, acp(3).*m1c, acp(4).*m1c];
    
    m2c = [0, av, v0-pi(2)];
    m2 = [acp(1).*m2c, acp(2).*m2c, acp(3).*m2c, acp(4).*m2c];
    
    M(2*i-1,:)=m1;
    M(2*i,:)=m2;
end

%slove MV = 0;
[V,D] = eig(M'*M);
D = max(D);
D(~D) = inf;
[~, i] = min(D);
V = V(:,i);

%solve the scale of control point
Cc = reshape(V,3,4);

vd12 = pdist2(Cc(:,1),Cc(:,2));
vd13 = pdist2(Cc(:,1),Cc(:,3));
vd14 = pdist2(Cc(:,1),Cc(:,4));
vd23 = pdist2(Cc(:,2),Cc(:,3));
vd24 = pdist2(Cc(:,2),Cc(:,4));
vd34 = pdist2(Cc(:,3),Cc(:,4));
vdall = [vd12 vd13 vd14 vd23 vd24 vd34];

cd12 = pdist2(Cw(:,1),Cw(:,2));
cd13 = pdist2(Cw(:,1),Cw(:,3));
cd14 = pdist2(Cw(:,1),Cw(:,4));
cd23 = pdist2(Cw(:,2),Cw(:,3));
cd24 = pdist2(Cw(:,2),Cw(:,4));
cd34 = pdist2(Cw(:,3),Cw(:,4));
cdall = [cd12 cd13 cd14 cd23 cd24 cd34];

scale = sum(vdall.*cdall) / sum(vdall.^2);
Cc = scale.*Cc;

%estimate the camera position
[R,T] = caculateRT(Cw,Cc);

end
