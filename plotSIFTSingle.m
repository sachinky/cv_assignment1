function plotSIFTSingle(img, points)
    pSize = size(points, 2);
    coords = zeros(pSize, 2);
    
    for i=1:pSize
        coords(i, 1) = points{i}.x;
        coords(i, 2) = points{i}.y;
    end
    
    figure;
    imshow(img);
    colormap(gray(256));
    axis image;
    hold on;
    plot(coords(:, 1), coords(:, 2), 'go');
end