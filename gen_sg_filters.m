clearvars; clc;

% generate Savitsky Golay filter coefficients for frame lengths of odd numbers from 1 to
% 300 or so and save in .mat file

max = 401;
Ns = (max-1)/2;

sg_smooth = cell(Ns,1);
sg_diff = cell(Ns,1);

for i = 1:Ns
    fl = (2*i) + 1;

    [b,g] = sgolay(2,fl);
    sg_smooth{i} = b;
    sg_diff{i} = g;
end


save('sg_filters.mat', 'sg_smooth', 'sg_diff');