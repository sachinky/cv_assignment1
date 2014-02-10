function [r, count] = SIFT(imPath)
    img = imread(imPath);
    gray  = im2double(rgb2gray(img));
    gss = gaussianScaleSpace(gray);
    DoG = dogScaleSpace(gss);
    sft = getSIFTFeatures(DoG);
    [points, c0] = getSIFTKeypoints(sft);
    [refinedPoints, c1] = refinePoints(DoG, points, c0);
    [orientedPoints, c2, gradM, gradN, normG, atanG] = getOrientedSIFTKeypoints(gss, refinedPoints, c1);
    [featureVectors, c3] = getSIFTFeatureVectors(gss, orientedPoints, c2, gradM, gradN, normG, atanG);
    r = featureVectors;
    count = c3;
%     r = orientedPoints;
end