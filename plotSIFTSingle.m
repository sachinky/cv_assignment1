function plotSIFTSingle(img, points, pSize)
%     pSize = size(points, 2);
    coords = zeros(pSize, 3);
    
%     for i=1:pSize
%         coords(i, 1) = points{i}.x;
%         coords(i, 2) = points{i}.y;
%         coords(i, 3) = 7.5*points{i}.sigma;
%     end
    
    figure;
    imshow(img);
%     colormap(gray(256));
    axis image;
    hold on;
%     plot(coords(:, 1), coords(:, 2), 'go');
    for i=1:pSize
        x = double(points{i}.x);
        y = double(points{i}.y);
        theta = points{i}.theta;
        r = 7.5*points{i}.sigma;
        circle(x, y, r);
        plot([x, x+r*cos(theta)], [y, y+r*sin(theta)], 'Color','r','LineWidth',2);
    end

end

function circle(x,y,r)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp);
end