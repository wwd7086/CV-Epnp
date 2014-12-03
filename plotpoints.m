function plotpoints(index,pts1,col2,surfaceNum,x3d,col3)

col2 = strcat(col2,'o');
col3 = strcat(col3,'o');
for i=1:numel(surfaceNum)
    hold on;
    subplot(1,5,index);
    plot(pts1(1,i),pts1(2,i),col2);
    switch surfaceNum(i)
        case 0
            return;
        case 1
            coord = [1,2];
        case 2
            coord = [3,2];
        case 3
            coord = [3,1];
    end
    subplot(1,5,surfaceNum(i)+1);
    axis([0 200 0 200]);
    axis ij;
    x3ds = x3d(:,i);
    hold on;
    plot(x3ds(coord(1)),x3ds(coord(2)),col3);
end

end