function w = dogScaleSpace(gss)
    n_oct = size(gss, 2);
    n_spo = 3;
    
    w = {};
    for i=1:n_oct
        w{i} = {};
        for j=1:n_spo+2
            w{i}{j} = gss{i}{j+1} - gss{i}{j};
        end
    end
end