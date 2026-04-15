function [b,r] = scale_analog(an,pc)
    % Basic linear regression
    % X = [ones(length(an),1) an];
    % b = X\pc; 
    % Correlation coefficient
    % r = corrcoef(an,pc);

    % Robust linear regression
    [b,stats] = robustfit(an,pc);
    r = stats.coeffcorr(1,2);

end