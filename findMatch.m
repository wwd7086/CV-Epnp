function [pts1, pts2] = findMatch(I1,I2,tresh,num)

%Get the key point
Options.upright=false;
Options.tresh=tresh; %0.0001
Ipts1=OpenSurf(I1,Options);
Ipts2=OpenSurf(I2,Options);

pts1 = zeros(2,num);
pts2 = zeros(2,num);

% Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]); 
% Find the best matches
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  for i=1:length(Ipts1)
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
      [err(i),cor2(i)]=min(distance);
  end
% Sort matches on vector distance
  [~, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
  
for i=1:num
    pts1(1,i)=Ipts1(cor1(i)).x;
    pts1(2,i)=Ipts1(cor1(i)).y;
    
    pts2(1,i) =Ipts2(cor2(i)).x;
    pts2(2,i) = Ipts2(cor2(i)).y;
end
end