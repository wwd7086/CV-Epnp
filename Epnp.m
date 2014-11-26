function [R,T] = Epnp(K,Pi,Pw)
%number of the control point
numCp = 3;

au = K(1,1);
av = K(2,2);
at = K(1,2);
u0 = K(1,3);
v0 = K(2,3);

numPw = size(Pw,2);

%pick 4 random points as control point
randCi = randperm(numPw);
randCi = randCi(1:numCp);
Cw = Pw(:,randCi);

%express all the point using control point
Acp = zeros(numCp,numPw);
for i=1:numPw

    acp = Cw\Pw(:,i);
    Acp(:,i) = acp;
    
end

%solve the control point in camer coordinate
M = zeros(numPw*2,3*numCp);
for i = 1:numPw
    acp = Acp(:,i);
    pi = Pi(:,i);
    
    m1c = [au, at, u0-pi(1)];
    m1 = kron(acp',m1c);
    %m1 = [acp(1).*m1c, acp(2).*m1c, acp(3).*m1c, acp(4).*m1c];
    
    m2c = [0, av, v0-pi(2)];
    m2 = kron(acp', m2c);
    %m2 = [acp(1).*m2c, acp(2).*m2c, acp(3).*m2c, acp(4).*m2c];
    
    M(2*i-1,:)=m1;
    M(2*i,:)=m2;
end

%slove MV = 0;
[s,d,V] = svd(M);
V = V(:,end);

%solve the scale of control point
Cc = reshape(V,3,numCp);
distCc = pdist2(Cc',Cc');
distCw = pdist2(Cw',Cw');
rescale = distCw./distCc;
rescale = rescale(rescale>0);
rescale = mean(mean(rescale));
Cc = rescale.*Cc;

%estimate the camera position
[R,T] = caculateRT(Cw,Cc);

end
