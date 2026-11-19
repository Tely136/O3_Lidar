function tests = test_resolution
    tests = functiontests(localfunctions);
end
function setupOnce(t)
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'proc_functions')));
end
function testIdentity(t)
    verifyEqual(t,resolution_IR(1,3.75),3.75);
    verifyEqual(t,resolution_FC(1,3.75),3.75);
end
function testBoxcar(t)
    verifyEqual(t,resolution_IR(ones(5,1),2),10,'AbsTol',1e-12);
    fc = fzero(@(f) sin(5*pi*f)/(5*sin(pi*f)) - 0.5,[0.01 0.19]);
    verifyEqual(t,resolution_FC(ones(5,1),2),1/fc,'RelTol',1e-6);
end
function testDerivative(t)
    c = [-0.5;0;0.5];
    verifyEqual(t,resolution_IR(c,2),4,'AbsTol',1e-12);
    fc = fzero(@(f) sin(2*pi*f)/(2*pi*f) - 0.5,[0.01 0.49]);
    verifyEqual(t,resolution_FC(c,2),1/fc,'RelTol',1e-6);
    verifyEqual(t,resolution_IR(-3*c,2),resolution_IR(c,2));
    verifyEqual(t,resolution_FC(-3*c,2),resolution_FC(c,2),'AbsTol',1e-12);
end
function testFilterCells(t)
    c = [ones(3,1)/3 [-0.5;0;0.5]];
    C = {c c};
    verifyEqual(t,resolution_IR(C,2),[4 4]);
    verifyEqual(t,resolution_FC(C,2),repmat(resolution_FC(c(:,2),2),1,2));
    verifyEqual(t,resolution_IR(c(:,1),2),6);
end
function testSavitzkyGolay(t)
    c = [-3;12;17;12;-3] / 35;
    verifyEqual(t,resolution_IR(c,1),37/15,'AbsTol',1e-12);
    fc = fzero(@(f) (17 + 24*cos(2*pi*f) - 6*cos(4*pi*f))/35 - 0.5,[0.1 0.4]);
    verifyEqual(t,resolution_FC(c,1),1/(2*fc),'RelTol',1e-6);
end
function testSuccessiveFilters(t)
    c = conv(ones(3,1)/3,ones(3,1)/3);
    verifyEqual(t,resolution_IR(c,1),3,'AbsTol',1e-12);
    fc = fzero(@(f) ((1 + 2*cos(2*pi*f))/3)^2 - 0.5,[0.01 0.3]);
    verifyEqual(t,resolution_FC(c,1),1/(2*fc),'RelTol',1e-6);
end
function testInvalidFilters(t)
    verifyError(t,@()resolution_IR([1 1],1),'O3:FilterLength');
    verifyError(t,@()resolution_FC([1 0 2],1),'O3:FilterSymmetry');
    verifyError(t,@()resolution_FC([0 0 0],1),'O3:FilterNormalization');
end
