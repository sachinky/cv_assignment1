function [r, count] = getSIFTFeatureVectors(gss, points, c2, normG, atanG)
    del_min = 0.5;
    n_ori = 8;
    n_hist = 4;
    lambda_descr = 6.0;
    
    c = 1;
    r = cell(c2);
    
    for it=1:c2
        fCount = 0;
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
            if m<1 || m>limM
                continue;
            end
            
            for n=minN:maxN
                if n<1 || n>limN
                    continue;
                end
                x_cap = ((m*delta_okey - x_key)*cos(theta_key) + (n*delta_okey - y_key)*sin(theta_key))/sigma_key;
                y_cap = (-(m*delta_okey - x_key)*sin(theta_key) + (n*delta_okey - y_key)*cos(theta_key))/sigma_key;
                
                if max(x_cap, y_cap)<lambda_descr*(n_hist+1)/n_hist
                    theta_cap = mod(atanG{o_key}{s_key}(n, m) - theta_key, 2*pi);
                    
                    c_descr = 1.0/(sqrt(2*pi)*lambda_descr*sigma_key)*exp(-double(((m*delta_okey - x_key)^2 + (n*delta_okey - y_key)^2)/(2*(lambda_descr*sigma_key)^2)))*normG{o_key}{s_key}(n, m);
                    
                    for i=1:n_hist
                        if abs(x_cent(i) - x_cap) <= 2*lambda_descr/n_hist
                            for j=1:n_hist
                                if abs(y_cent(j) - y_cap) <= 2*lambda_descr/n_hist
                                    for k=1:n_ori
                                        if mod(theta_cap - (2*pi)*(k-1)/n_ori, 2*pi) < 2*pi/n_ori
                                            fCount = fCount + 1;
                                            h(i, j, k) = h(i, j, k) + (1 - n_hist/(2*lambda_descr)*abs(x_cap - x_cent(i)))*(1 - n_hist/(2*lambda_descr)*abs(y_cap - y_cent(j)))*(1 - n_ori/(2*pi)*abs(mod(theta_cap - theta_key - (2*pi)*(k-1)/n_ori, 2*pi)))*1000*c_descr;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        fCount
        f = zeros(1, n_hist*n_hist*n_ori);
        
        for i=1:n_hist
            for j=1:n_hist
                for k=1:n_ori
                    f((i-1)*n_hist*n_ori + (j-1)*n_ori + k) = h(i, j, k);
                end
            end
        end
        
        fNorm = norm(f);
        f = min(f, 0.2*fNorm);
        
        fNorm = norm(f);
        f = f./fNorm;
        
        f = uint8(256.*f);
        
        r{c} = struct('y', y_key, 'x', x_key, 'sigma', sigma_key, 'octave', o_key, 's', s_key, 'theta', theta_key, 'feature', f);
        c = c + 1;
    end
    count = c - 1;
end