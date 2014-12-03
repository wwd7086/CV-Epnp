function [R,T] = Epnp(K,Pi,Pw)
%number of the control point
numCp = 4;

au = K(1,1);
av = K(2,2);
at = K(1,2);
u0 = K(1,3);
v0 = K(2,3);

numPw = size(Pw,2);

isValid = false;
Acp = zeros(numCp,numPw);
while ~isValid 
    %pick 4 random points as control point
    randCi = randperm(numPw);
    randCi = randCi(1:numCp);
    Cw = Pw(:,randCi);
    
    if numCp>=4
        Cwt = [Cw;[1,1,1,1]];
        %express all the point using control point
        for i=1:numPw
            acp = Cwt\[Pw(:,i);1];
            Acp(:,i) = acp;  
        end
    elseif numCp==3
        %express all the point using control point
        for i=1:numPw
            acp = Cw\Pw(:,i);
            Acp(:,i) = acp;  
        end
    end
    
    if isempty(find(Acp==Inf)) && ~any(any(isnan(Acp)))     
        if numCp>=4 && rank(Cw)==3
            isValid = true;

        elseif numCp==3 && rank(Cw)>=2
            isValid = true;
        end
    end
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
distCc = sqrt(pdist2(Cc',Cc'));
distCw = sqrt(pdist2(Cw',Cw'));
rescale = distCw./distCc;

%eliminate abnormal value
rescale = rescale(rescale>0);
rescale(rescale == Inf) = 10000000;

rescale = mean(mean(rescale));
Cc = rescale.*Cc;

%estimate the camera position
[R,T] = caculateRT(Cw,Cc);

end
