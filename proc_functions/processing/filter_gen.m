function C = filter_gen(M1,M2,h1,h2,hkm)

    % Determine filter widths and load filter coefficients
    fl = gen_framelength(M1,M2,h1,h2,hkm);
    
    % TODO: don't load filters here, make altitude dependent filter list
    % previously and input to processsing function
    sg = load('sg_filters.mat');
    sg_diff = sg.sg_diff;

    C = cell(1,length(fl));
    for k = 1:length(fl)
        m = fl(k);
        N = (m-1)/2;

        C{k} = sg_diff{N};
    end
end