function y = deadtime_corr(x,deadt)
y = x./(1-x./deadt);
end