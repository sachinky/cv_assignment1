function plotSIFT(img1, img2, matching)
    mSize = size(matching, 2);
    points1 = zeros(mSize, 2);
    points2 = zeros(mSize, 2);
    
    for i=1:mSize
        points1(i, :) = [matching{i}(1).x matching{i}(1).y];
        points2(i, :) = [matching{i}(2).x matching{i}(2).y];
    end
    
    figure; showMatchedFeatures(rgb2gray(imread(img1)), rgb2gray(imread(img2)), points1, points2);
end