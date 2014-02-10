function r = getSIFTFeatures(DoG)
    n_spo = 3;
    n_oct = 4;
    c_DoG = 0.03;
    c_edge = 10.0;
    edgeThres = (c_edge+1)*(c_edge+1)/c_edge;
    
    r = cell(n_oct, n_spo);
    for oct = 1:n_oct
%         r{oct} = {};
        height = size(DoG{oct}{1}, 1);
        width = size(DoG{oct}{1}, 2);
        for spo = 2:n_spo+1
            A = DoG{oct}{spo} > imdilate(DoG{oct}{spo}, [1 1 1; 1 0 1; 1 1 1]);
            B = imdilate(DoG{oct}{spo-1}, [1 1 1; 1 1 1; 1 1 1]);
            C = imdilate(DoG{oct}{spo+1}, [1 1 1; 1 1 1; 1 1 1]);
            D = imdilate(DoG{oct}{spo}, [1 1 1; 1 0 1; 1 1 1]);
            
            Am = -DoG{oct}{spo} > imdilate(-DoG{oct}{spo}, [1 1 1; 1 0 1; 1 1 1]);
            Bm = imdilate(-DoG{oct}{spo-1}, [1 1 1; 1 1 1; 1 1 1]);
            Cm = imdilate(-DoG{oct}{spo+1}, [1 1 1; 1 1 1; 1 1 1]);
            Dm = imdilate(-DoG{oct}{spo}, [1 1 1; 1 0 1; 1 1 1]);
            
            H = zeros(height, width, 3); %H is the 2-D Hessian matrix
            H(:, :, 1) = imfilter(DoG{oct}{spo}, [1 -2 1]);
            H(:, :, 2) = imfilter(DoG{oct}{spo}, [0.25 0 -0.25; 0 0 0; -0.25 0 0.25]);
            H(:, :, 3) = imfilter(DoG{oct}{spo}, [1; -2; 1]);
            
            edgeness = ((H(:, :, 1) + H(:, :, 3)) .* (H(:, :, 1) + H(:, :, 3))) ./ (H(:, :, 1).*H(:, :, 3) - H(:, :, 2).*H(:, :, 2));
%             edgeness = zeros(height, width);
            r{oct}{spo-1} = abs(DoG{oct}{spo}) > 0.8*c_DoG & edgeness < edgeThres & ((A & D > max(B, C)) | (Am & Dm > max(Bm, Cm)));
        end
    end
    
end