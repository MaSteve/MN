function [sol, flag] = jacobi(A, b, iter, prec)
    flag = 0;
    n = size(A)(1);
    normb = norm(b);
    r = zeros(n, 1);
    d = zeros(n, 1);
    u_k = zeros(n, 1);
    for k = 1:iter
        r = b' - A*u_k;
        if prec*normb > norm(r)
            flag = 1;
            break;
        end
        d = r./diag(A);
        u_k = u_k + d;
    end
    sol = u_k;
end