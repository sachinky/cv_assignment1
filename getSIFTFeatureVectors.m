function [r, count] = getSIFTFeatureVectors(gss, points, c2)
    n_oct = 4;
    n_spo = 3;
    del_min = 0.5;
    n_ori = 8;
    n_hist = 4;
    lambda_descr = 6.0;
    
    gradM = gradient(gss, n_oct, n_spo, 'x');
    gradN = gradient(gss, n_oct, n_spo, 'y');
    
    c = 1;
    r = cell(c2);
    
    for it=1:c2
%         count = 0;
        h = zeros(n_hist, n_hist, n_ori);
        x_key = points{it}.x;
        y_key = points{it}.y;
        sigma_key = points{it}.sigma;
        o_key = points{it}.octave;
        s_key = points{it}.s;
        delta_okey = del_min*realpow(2, points{it}.octave-1);
        theta_key = points{it}.theta;
        
        maxM = uint16(round(x_key + sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/delta_okey);
        minM = uint16(round(x_key - sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/delta_okey);
        maxN = uint16(round(y_key + sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/delta_okey);
        minN = uint16(round(y_key - sqrt(2)*lambda_descr*sigma_key*(n_hist+1)/n_hist)/delta_okey);
        
        limM = size(gss{o_key}{1}, 2);
        limN = size(gss{o_key}{1}, 1);
        
        x_cent = zeros(1, n_hist);
        y_cent = zeros(1, n_hist);
        
        for i=1:n_hist
            x_cent(i) = (i - (1+n_hist)/2.0)*2*lambda_descr/n_hist;
            y_cent(i) = (i - (1+n_hist)/2.0)*2*lambda_descr/n_hist;
        end
        
        for m=minM:maxM
            for n=minN:maxN
                x_cap = ((m*delta_okey - x_key)*cos(theta_key) + (n*delta_okey - y_key)*sin(theta_key))/sigma_key;
                y_cap = (-(m*delta_okey - x_key)*sin(theta_key) + (n*delta_okey - y_key)*cos(theta_key))/sigma_key;
                
                if (max(x_cap, y_cap)<lambda_descr*(n_hist+1)/n_hist && m>=1 && n>=1 && m<=limM && n<=limN)
                    theta_cap = mod(atan2(gradM{o_key}{s_key}(n, m), gradN{o_key}{s_key}(n, m)) - theta_key, 2*pi);
                    
                    c_descr = 1.0/(sqrt(2*pi)*lambda_descr*sigma_key)*exp(-double(((m*delta_okey - x_key)^2 + (n*delta_okey - y_key)^2)/(2*(lambda_descr*sigma_key)^2)))*((gradM{o_key}{s_key}(n, m) - gradN{o_key}{s_key}(n, m))^2);
                    
                    for i=1:n_hist
                        for j=1:n_hist
                            if (abs(x_cent(i) - x_cap) <= 2*lambda_descr/n_hist && abs(x_cent(j) - y_cap) <= 2*lambda_descr/n_hist)
                                for k=1:n_ori
                                    if mod(theta_cap - theta_key - (2*pi)*(k-1)/n_ori, 2*pi) < 2*pi/n_ori
                                        
%                                         if c_descr>0.01
%                                             count = count+1;
%                                         end
                                        
                                        h(i, j, k) = h(i, j, k) + (1 - n_hist/(2*lambda_descr)*abs(x_cap - x_cent(i)))*(1 - n_hist/(2*lambda_descr)*abs(y_cap - y_cent(j)))*(1 - n_ori/(2*pi)*abs(mod(theta_cap - theta_key - (2*pi)*(k-1)/n_ori, 2*pi)))*1000*c_descr;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
%         count
        f = zeros(1, (n_hist^2)*n_ori);
        
        for i=1:n_hist
            for j=1:n_hist
                for k=1:n_ori
                    f((i-1)*n_hist*n_ori + (j-1)*n_ori + k) = h(i, j, k);
                end
            end
        end
        
        fNorm = norm(f);
        
        f = min(f./fNorm, 0.2);
        f = uint8(256.*f);
        
        r{c} = struct('y', y_key, 'x', x_key, 'sigma', sigma_key, 'octave', o_key, 's', s_key, 'theta', theta_key, 'feature', f);
        c = c + 1;
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