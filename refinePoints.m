function [r, count] = refinePoints(DoG, points, c0)
    n_spo = 3;
    delMin = 0.5;
    sigmaMin = 0.8;
    c_DoG = 0.03;

    r = cell(c0);
    c = 1;
    for i = 1:c0
        m = points{i}.m;
        n = points{i}.n;
        sigma_key = points{i}.sigma;
        o_key = points{i}.octave;
        s_key = points{i}.s;
        delta_okey = delMin*realpow(2, points{i}.octave-1);
        
        limM = size(DoG{o_key}{1}, 2);
        limN = size(DoG{o_key}{1}, 1);
        
        for iter = 1:5
            temp_mat = zeros(3,3,3);
            for a = -1:1
                for b = -1:1
                    for d = -1:1
                        temp_mat(2+a,2+b,2+d) = DoG{o_key}{s_key+a}(n+d, m+b);
                    end
                end
            end
        
            grad = zeros(3,1);
            H = zeros(3,3);
        
            grad(1,1) = (temp_mat(3,2,2) - temp_mat(1,2,2))/2;
            grad(2,1) = (temp_mat(2,3,2) - temp_mat(2,1,2))/2;  
            grad(3,1) = (temp_mat(2,2,3) - temp_mat(2,2,1))/2;

            H(1,1) = temp_mat(3,2,2) + temp_mat(1,2,2) - 2*temp_mat(2,2,2);  
            H(2,2) = temp_mat(2,3,2) + temp_mat(2,1,2) - 2*temp_mat(2,2,2);  
            H(3,3) = temp_mat(2,2,3) + temp_mat(2,2,1) - 2*temp_mat(2,2,2);  
            H(1,2) = (temp_mat(3,3,2) - temp_mat(3,1,2) - temp_mat(1,3,2) + temp_mat(1,1,2))/4;
            H(2,1) = H(1,2);
            H(1,3) = (temp_mat(3,2,3) - temp_mat(3,2,1) - temp_mat(1,2,3) + temp_mat(1,2,1))/4;
            H(3,1) = H(1,3); 
            H(2,3) = (temp_mat(2,3,3) - temp_mat(2,3,1) - temp_mat(2,1,3) + temp_mat(2,1,1))/4;
            H(3,2) = H(2,3);

            H_inv = inv(H);
            alpha = - H_inv * grad;
            if(max(abs(alpha)) < 0.6)
                omega = temp_mat(2,2,2) - 0.5* grad' * H_inv *grad;
                if(abs(omega) > c_DoG)
                    sigma = s_key*sigmaMin*realpow(2, (sigma_key+alpha(1))/n_spo);
                    x = round(delta_okey * (m + alpha(2)));
                    y = round(delta_okey * (n + alpha(3)));
                    r{c} = struct('y', y, 'x', x, 'sigma', sigma, 'octave', o_key, 's', s_key);
                    c = c + 1;
                    break;
                end
            else
                s_key = min(max(2, uint8(round(s_key + alpha(1)))), n_spo+1);
                m = min(max(2, uint16(round(m + alpha(2)))), limM-1);
                n = min(max(2, uint16(round(n + alpha(3)))), limN-1);
            end
        end
    end
    count = c - 1;
end