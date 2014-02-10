function r = SIFT(imPath)
    img = imread(imPath);
    gray  = im2double(rgb2gray(img));
    gss = gaussianScaleSpace(gray);
    DoG = dogScaleSpace(gss);
    sft = getSIFTFeatures(DoG);
    points = getSIFTKeypoints(sft);
    orientedPoints = getOrientedSIFTKeypoints(gss, points);
    featureVectors = getSIFTFeatureVectors(gss, orientedPoints);
    r = featureVectors;
%     r = orientedPoints;
end