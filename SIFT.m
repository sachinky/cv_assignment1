function [r, count] = SIFT(imPath)
    img = imread(imPath);
    gray  = im2double(rgb2gray(img));
    gss = gaussianScaleSpace(gray);
    DoG = dogScaleSpace(gss);
    sft = getSIFTFeatures(DoG);
    [points, c0] = getSIFTKeypoints(sft);
    c0
    [refinedPoints, c1] = refinePoints(DoG, points, c0);
    c1
    [orientedPoints, c2, normG, atanG] = getOrientedSIFTKeypoints(gss, refinedPoints, c1);
    c2
    [featureVectors, c3] = getSIFTFeatureVectors(gss, orientedPoints, c2, normG, atanG);
    c3
    r = featureVectors;
    count = c3;
%     r = orientedPoints;
end