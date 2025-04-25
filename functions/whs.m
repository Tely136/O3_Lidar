function z = whs(y,lambda)
m = length(y);
if m<1000
    E = eye(m);
    D = diff(E);
    z = (E + lambda * D' * D)\y;
else
    % Invoke sparse method for large m
    E = speye(m);
    D = diff(E);
    C = chol(E + lambda * D' * D);
    z = C\(C'\y);
end