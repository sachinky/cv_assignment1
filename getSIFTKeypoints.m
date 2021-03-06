function [r, count] = getSIFTKeypoints(pss)
    n_oct = 4;
    n_spo = 3;
    delMin = 0.5;
    sigmaMin = 0.8;
    
    total = 0;
    for i=1:n_oct
        for j=1:n_spo
            total = total + sum(sum(pss{i}{j}));
        end
    end
    r = cell(total);
    c = 1;
    
    imHeight = size(pss{1}{1}, 1)*delMin;
    imWidth = size(pss{1}{1}, 2)*delMin;
    

    for oct=1:n_oct
        height = size(pss{oct}{1}, 1);
        width = size(pss{oct}{1}, 2);
        for spo=1:n_spo
            for i=1:height
                for j=1:width
                    if pss{oct}{spo}(i, j)
                        scale = realpow(2, oct-1);
                        delta = delMin*scale;
                        sigma = scale*sigmaMin*realpow(2, (spo+1)/n_spo);
                        if not(uint16(i*delta) < 0.01*imHeight || uint16(i*delta) > 0.99*imHeight || uint16(j*delta) < 0.01*imWidth || uint16(j*delta) > 0.99*imWidth)
                            r{c} = struct('m', j, 'n', i, 'sigma', sigma, 'octave', oct, 's', spo+1);
                            c = c + 1;
                        end
                    end
                end
            end
        end
    end
    count = c - 1;
end