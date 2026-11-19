function tests = test_aerosol_correction
    tests = functiontests(localfunctions);
end

function testSyntheticReturns(t)
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'proc_functions')));

    dz = 15;
    z = (1200:dz:6000)';
    beta_m = 3e-5 .* exp(-z/8000);
    beta_a = 2e-5 .* exp(-((z-2600)/400).^2);
    beta_a(end) = 0;
    ozone = 1.5e18*ones(size(z));
    sigma_off = 5e-23*ones(size(z));
    dsigma = 2e-22*ones(size(z));
    lam_on = 287.2e-9;
    lam_off = 299.1e-9;
    molecular_scale = (lam_off/lam_on)^4;
    aerosol_scale = (lam_off/lam_on)^0.5;
    alpha_m = beta_m.*(8*pi/3);
    alpha_a = 60*beta_a;
    P_off = (beta_m + beta_a) ./ z.^2 .* exp(-2*cumtrapz(z,alpha_m + alpha_a + sigma_off.*ozone));
    P_on = (molecular_scale*beta_m + aerosol_scale*beta_a) ./ z.^2 .* ...
        exp(-2*cumtrapz(z,molecular_scale*alpha_m + aerosol_scale*alpha_a + (sigma_off+dsigma).*ozone));

    C = repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,numel(z));
    baseline = dial(P_on,P_off,C,dz,dsigma) - (molecular_scale-1)*alpha_m./dsigma;

    ids = 2:numel(z)-1;
    settings = struct('start_km',z(ids(1))/1e3,'reference_km',z(ids(end))/1e3,...
        'reference_beta',beta_a(ids(end)));
    [D_fr,D_nr,retrieved] = aerosol_correction(P_off,baseline,z,C,C,alpha_m,sigma_off,...
        dsigma,lam_on,lam_off,settings);

    corrected = baseline + D_fr;
    usable = ids(2:end-1);
    verifyTrue(t,all(isfinite(corrected(usable))));
    verifyLessThan(t,norm(corrected(usable)-ozone(usable)),0.02*norm(baseline(usable)-ozone(usable)));
    verifyLessThan(t,norm(retrieved(ids)-beta_a(ids)),0.02*norm(beta_a(ids)));
    verifyTrue(t,isequaln(D_fr,D_nr));
    verifyTrue(t,all(isnan(D_fr([1:2 end-1:end]))));
end
