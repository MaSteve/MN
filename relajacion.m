function [sol, flag] = relajacion(A, b, iter, prec, w)
    flag = 0;
    n = size(A)(1);
    normb = norm(b);
    r = zeros(n, 1);
    d = zeros(n, 1);
    u_k = zeros(n, 1);
    u_k1 = zeros(n, 1);
    for k = 1:iter
        r(1) = b(1) - A(1, 1:n)*u_k(1:n);
        d(1) = w*(r(1)/A(1, 1));
        u_k1(1) = u_k(1) + d(1);
        for i = 2: n
            r(i) = b(i) - A(i,1:i-1)*u_k1(1:i-1) - A(i, i:n)*u_k(i:n);
            d(i) = w*(r(i)/A(i, i));
            u_k1(i) = u_k(i) + d(i);
        end
        if prec*normb > norm(r)
            flag = 1;
            break;
        end
        u_k = u_k1;
    end
    sol = u_k1;
end
