function [r, count] = refinePoints(DoG, points, c0)
    n_spo = 3;
    delMin = 0.5;
    sigmaMin = 0.8;
    c_DoG = 0.03;

    c0
    r = cell(c0);
    c = 1;
    for i = 1:c0
%         x_key = points{i}.x;
%         y_key = points{i}.y;
        m = points{i}.m;
        n = points{i}.n;
        sigma_key = points{i}.sigma;
        o_key = points{i}.octave;
        s_key = points{i}.s;
        delta_okey = delMin*realpow(2, points{i}.octave-1);
        
        for iter = 1:1
            temp_mat = zeros(3,3,3);
            for a = -1:1
                for b = -1:1
                    for c = -1:1
                        temp_mat(2+a,2+b,2+c) = DoG{o_key}{s_key}(n+c, m+b);
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
            if(max(abs(alpha)) <0.5)
                omega = temp_mat(2,2,2) - 0.5* grad' * H_inv *grad;
                if(abs(omega) > c_DoG)
%                     delta_o = 2^(oct-1) / delta_min;
%                     new_s =  (delta_o*sigma_min/delta_min)*(2^(alpha(1)+scale)/n_spo);
                    sigma = s_key*sigmaMin*realpow(2, (sigma_key+alpha(1))/n_spo);
%                     new_m = delta_o *(alpha(2) + m);
%                     new_n = delta_o *(alpha(3) + n);
%                     point = [oct scale m n new_m new_n new_s omega];
%                     refined = [refined;point];
                    x = uint16(delta_okey * (m + alpha(2)));
                    y = uint16(delta_okey * (n + alpha(3)));
                    r{c} = struct('y', y, 'x', x, 'sigma', sigma, 'octave', o_key, 's', s_key);
                    c = c + 1;
                    break;
                end
%             else
%                 scale = round(scale + alpha(1));
%                 m = round(m + alpha(2));
%                 n = round(n + alpha(3));
            end
        end
    end
    count = c - 1;
end