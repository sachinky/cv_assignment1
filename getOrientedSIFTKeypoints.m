function [r, count] = getOrientedSIFTKeypoints(gss, points, c1)
    lambda_ori = 1.5;
    n_bins = 36;
    n_oct = 4;
    n_spo = 3;
    del_min = 0.5;
    thres_hist = 0.8;
    
    gradM = gradient(gss, n_oct, n_spo, 'x');
    gradN = gradient(gss, n_oct, n_spo, 'y');
    
    c=1;
    
    r = cell(2*c1);
    
    for i=1:c1
        x_key = points{i}.x;
        y_key = points{i}.y;
        sigma_key = points{i}.sigma;
        o_key = points{i}.octave;
        s_key = points{i}.s;
        delta_okey = del_min*realpow(2, points{i}.octave-1);
        maxM = uint16(round(x_key + 3*lambda_ori*sigma_key)/delta_okey);
        minM = uint16(round(x_key - 3*lambda_ori*sigma_key)/delta_okey);
        maxN = uint16(round(y_key + 3*lambda_ori*sigma_key)/delta_okey);
        minN = uint16(round(y_key - 3*lambda_ori*sigma_key)/delta_okey);
%         
%         if(minM<1 || minN<1 || maxM>size(gss{o_key}{1}, 2) || maxN>size(gss{o_key}{1}, 1))
%             continue;
%         end
        
        h_cur = zeros(1, n_bins);
        limM = size(gss{o_key}{1}, 2);
        limN = size(gss{o_key}{1}, 1);
        
        for m=minM:maxM
            
            if m<1 || m>limM
                continue
            end
            
            for n=minN:maxN
                
                if n<1 || n>limN
                    continue
                end
                
                c_ori = 1/(sqrt(2*pi)*lambda_ori*sigma_key)*exp(-double(((m*delta_okey - x_key)^2 + (n*delta_okey - y_key)^2)/(2*(lambda_ori*sigma_key)^2)))*((gradM{o_key}{s_key}(n, m) - gradN{o_key}{s_key}(n, m))^2);
                index = uint16(round(n_bins/(2*pi)*atan2(gradM{o_key}{s_key}(n, m), gradN{o_key}{s_key}(n, m)))) + 1;
                
                h_cur(index) = h_cur(index) + c_ori;
            end
        end
        
        h_new = zeros(1, n_bins);
        for j=1:6
            for bin=1:n_bins
                if bin==1
                    h_new(bin) = h_cur(n_bins) + h_cur(1) + h_cur(2);
                else
                    if bin==n_bins
                        h_new(bin) = h_cur(n_bins-1) + h_cur(n_bins) + h_cur(1);
                    else
                         h_new(bin) = h_cur(bin-1) + h_cur(bin) + h_cur(bin+1);
                    end
                end
                
                h_new(bin) = h_cur(bin)/3.0;
            end
            
            h_cur = h_new;
        end
        
        for bin=1:n_bins
            h_km = h_cur(mod(bin-2, n_bins)+1);
            h_kM = h_cur(mod(bin, n_bins)+1);
            if h_cur(bin)>h_km && h_cur(bin)>h_kM && h_cur(bin) > thres_hist*max(h_cur)
                theta_key = (2*pi)*(bin-1)/n_bins + pi/n_bins*((h_km - h_kM)/(h_km + 2*h_cur(bin) + h_kM));
                r{c} = struct('y', y_key, 'x', x_key, 'sigma', sigma_key, 'octave', o_key, 's', s_key, 'theta', theta_key);
                c = c + 1;
            end
        end
    end
    count = c - 1;
end


function g = gradient(dog, n_oct, n_spo, dir)
    g = cell(n_oct, n_spo);
    for i=1:n_oct
%         g{i} = {};
        for j=1:n_spo+3
            if dir=='x'
                g{i}{j} = imfilter(dog{i}{j}, [-0.5 0 0.5]);
            else
                g{i}{j} = imfilter(dog{i}{j}, [-0.5; 0; 0.5]);
            end
        end
    end
end