function v = gaussianScaleSpace(imgIn)
    delMin = 0.5;
    sigmaIn = 0.5;
    sigmaMin = 0.8;
    n_spo = 3;
    n_oct = log2(min(size(imgIn, 1), size(imgIn, 2))) - 4;
    
    initImg = imresize(imgIn, 1/delMin, 'bilinear');
    v = {};
    
    sigma = 1/delMin * sqrt(sigmaMin*sigmaMin - sigmaIn*sigmaIn);
    curFilter = getGaussian(sigma);
    
    v{1} = {};
    v{1}{1} = imfilter(initImg, curFilter);
    
    for i=2:n_spo+3
        sigma = sigmaMin/delMin * sqrt(realpow(2, 2*(i-1)/n_spo) - realpow(2, 2*(i-2)/n_spo));
        curFilter = getGaussian(sigma);
        v{1}{i} = imfilter(v{1}{i-1}, curFilter);
    end
    
    for j=2:n_oct
        v{j} = {};
        v{j}{1} = imresize(v{j-1}{1}, 0.5);
        
        for i=2:n_spo+3
            sigma = sigmaMin/delMin * sqrt(realpow(2, 2*(i-1)/n_spo) - realpow(2, 2*(i-2)/n_spo));
            curFilter = getGaussian(sigma);
            v{j}{i} = imfilter(v{j}{i-1}, curFilter);
        end
    end
end

function r = getGaussian(sigma)
    r = fspecial('gaussian', floor(4*sigma), sigma);
end