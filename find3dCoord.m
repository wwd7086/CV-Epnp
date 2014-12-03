function [x3d, surfaces] = find3dCoord(masks, H, x2d, objectSize)
% masks is a 1*3 cell array containing
numPoint = size(x2d,2);
x3d = zeros(3,numPoint);
surfaces = zeros(1,numPoint);

x2d = round(x2d);
for i=1:numPoint
    
    %---- find belonged surface -----
    tag(1) = masks{1}(x2d(2,i),x2d(1,i))>0;
    tag(2) = masks{2}(x2d(2,i),x2d(1,i))>0;
    tag(3) = masks{3}(x2d(2,i),x2d(1,i))>0;

    surface = 0; 
    if tag(1)==1
        surface = 1;
    else if tag(2)==1
            surface = 2;
        else if tag(3)==1;
                surface = 3;
            end
        end
    end

    %---- find coordinates -----------
    switch surface
        case 0 % Not a feature point on the target object!
        x3ds = zeros(3,1); % Return empty
        case 1 % On the FRONT surface
            x2dhomo = [x2d(:,i)',1]';
            x3d_planar = H{1}*x2dhomo;
            x3d_planar = x3d_planar(1:2)*1/x3d_planar(3);
            x3ds(1:2) = bsxfun(@times,x3d_planar,objectSize(1:2)'/200);
            x3ds(3) = 0;
        case 2 % On the SIDE surface
            x2dhomo = [x2d(:,i)',1]';
            x3d_planar = H{2}*x2dhomo;
            x3d_planar = x3d_planar(1:2)*1/x3d_planar(3);
            x3ds([3,2]) = bsxfun(@times,x3d_planar,objectSize([3,2])'/200);
            x3ds(1) = objectSize(1);
        case 3 % On the TOP surface
            x2dhomo = [x2d(:,i)',1]';
            x3d_planar = H{3}*x2dhomo;
            x3d_planar = x3d_planar(1:2)*1/x3d_planar(3);
            x3ds([3,1]) = bsxfun(@times,x3d_planar,objectSize([3,1])'/200);
            x3ds(2) = 0;
    end

    if ~isempty(x3ds) && surface~=0
        x3ds = bsxfun(@minus, x3ds, objectSize/2);
        x3ds = x3ds';
    end
    
    x3d(:,i) = x3ds;
    surfaces(i) = surface;
    x3ds = [];

    
end



