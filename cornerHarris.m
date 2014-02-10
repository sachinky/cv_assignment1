function r=cornerHarris(imgPath, wSize, gaussian, sigma)
    img = imread(imgPath);
    img = rgb2gray(im2double(img));
    height = size(img, 1);
    width = size(img, 2);
    threshold = 0.1;
    
    gradientX = sobel(img, 'x');
    gradientY = sobel(img, 'y');
    
    H = zeros(height, width, 3);
    % For every pixel p, H(p) = [Gx(p)*Gx(p), Gx(p)*Gy(p); Gy(p)*Gx(p), Gy(p)*Gy(p)]
    % But since 2nd and 3rd elements of matrix are the same, we store them
    % only once
    H(:, :, 1) = gradientX .* gradientX;
    H(:, :, 2) = gradientX .* gradientY;
    H(:, :, 3) = gradientY .* gradientY;
    
%     if wSize>3
%         wSize = 3;
%     end
    if gaussian == false
        window = double(ones(2*wSize+1, 2*wSize+1));
    else
        window = fspecial('gaussian', [2*wSize+1, 2*wSize+1], sigma);
    end
    
    A = zeros(size(H));
    A = imfilter(H, window); %A is convolution of window with H matrix
    
    R = zeros(height, width); %R is the response matrix: R(p) = det(A(p)/trace(A(p))
    for i=1:height
        for j=1:width
            R(i, j) = (A(i, j, 1)*A(i, j, 3) - A(i, j, 2)*A(i, j, 2))/(A(i, j, 1) + A(i, j, 3) + 0.0000001);
            %R(i, j) = (A(i, j, 1)*A(i, j, 3) - A(i, j, 2)*A(i, j, 2)) - 0.05*(A(i, j, 1) + A(i, j, 3))*(A(i, j, 1) + A(i, j, 3));
        end
    end
    
    R = (R-min(min(R)))/(max(max(R))-min(min(R)));
    r1 = R > imdilate(R, [1 1 1; 1 0 1; 1 1 1]);
    r = r1 & (R > threshold);

end

function r=sobel(img, direction)
    if direction=='x'
        kernel = [1, 0, -1; 2, 0, -2; 1, 0, -1];
    else
        kernel = [1, 2, 1; 0, 0, 0; -1, -2, -1];
    end
    r=imfilter(img, kernel);
end