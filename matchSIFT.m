function [m, count] = matchSIFT(fA, fB)
    numFeaturesA = size(fA, 2);
    numFeaturesB = size(fB, 2);
    c_match_rel = 0.9;
    
    c = 1;
    m = cell(min(numFeaturesA, numFeaturesB));
    used = zeros(numFeaturesB, 1);
    
    for i=1:numFeaturesA
        small = 128*256;
        small2 = 128*256;

        index = 0;
        
%         if sqrt(sum(double(fA{i}.feature)))<20 || sqrt(sum(double(fA{i}.feature)))>51*120
%             continue;
%         end
        
        for j=1:numFeaturesB
            d = sqrt(sum((double(fA{i}.feature) - double(fB{j}.feature)).^2));
            if d<small
                small = d;
                index = j;
            else
                if d<small2
                    small2 = d;
                end
            end
        end
        
        if(small < 50 && small <= c_match_rel*small2 && used(index)==0)
            m{c} = [fA{i} fB{index}];
            used(index) = 1;
            c = c + 1;
        end
    end
    count = c - 1;
end