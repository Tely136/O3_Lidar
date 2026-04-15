function [b,r] = scale_analog(an,pc)
    % Basic linear regression
    % X = [ones(length(an),1) an];
    % b = X\pc; 
    % Correlation coefficient
    % r = corrcoef(an,pc);

    % Robust linear regression
    mdl = fitlm(an,pc,'linear','RobustOpts','on');
    % mdl = fitlm(an,pc,'linear');

    b = mdl.Coefficients.Estimate;
    r = sqrt(mdl.Rsquared.Ordinary);
end