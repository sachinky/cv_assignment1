function [r, count] = SIFT(imPath)
    img = imread(imPath);
    gray  = im2double(rgb2gray(img));
    gss = gaussianScaleSpace(gray);
    DoG = dogScaleSpace(gss);
    sft = getSIFTFeatures(DoG);
    [points, c1] = getSIFTKeypoints(sft);
    [orientedPoints, c2] = getOrientedSIFTKeypoints(gss, points, c1);
    [featureVectors, c3] = getSIFTFeatureVectors(gss, orientedPoints, c2);
    r = featureVectors;
    count = c3;
%     r = orientedPoints;
end