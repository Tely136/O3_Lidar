function tests = test_cloud_screening
    tests = functiontests(localfunctions);
end

function testCloudLayers(t)
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'proc_functions')));

    dz = 3.75;
    z = (dz:dz:6000)';
    clear_return = exp(-0.0005*z)./z.^2;
    enhancement = 3*exp(-((z-1100)/60).^2) + 4*exp(-((z-4500)/80).^2);
    far = clear_return .* [exp(enhancement) ones(size(z))];
    near = 4*far;
    overlap = find(z >= 1500 & z <= 2000);
    settings = struct('range_km',[0.5 5.5],'frame_length',11,...
        'positive_threshold',0.01,'negative_threshold',0.01);

    cloud = cloud_screening(near,far,z,overlap([1 end]),settings);
    layers = cloud.layers_km{1};
    verifySize(t,layers,[2 2]);
    verifyTrue(t,all(layers(:,1) < [1.1;4.5] & layers(:,2) > [1.1;4.5]));
    verifyFalse(t,any(cloud.mask(:,2)));

    C = repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,numel(z));
    ozone = dial(far,far,C,dz,1,cloud.mask);
    base = find(cloud.mask(:,1),1);
    verifyTrue(t,isnan(ozone(base-1,1)));
    verifyEqual(t,ozone(base-2,1),0);
end
