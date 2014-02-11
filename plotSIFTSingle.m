function plotSIFTSingle(img, points, pSize)
    sigmaIn = 0.5;
    
    figure;
    imshow(img);
    axis image;
    hold on;
    
    for i=1:pSize
        x = double(points{i}.x);
        y = double(points{i}.y);
        theta = points{i}.theta;
        r = 6.0*points{i}.sigma*sigmaIn;
        circle(x, y, r);
        text(max(1, x-2), y, '\fontsize{12}\color{blue}+');
        plot([x, x+r*cos(theta)], [y, y+r*sin(theta)], 'Color','b','LineWidth',1);
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