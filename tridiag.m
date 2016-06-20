function sol = tridiag(A, b) %Resuelve sistemas tridiagonales.
    n = size(A)(1);

    m = zeros(1,n);
    g = zeros(1,n);

    m(1) = A(1,1);
    g(1) = b(1)/m(1);

    for k = 2:n
        m(k) = A(k,k) - (A(k-1,k)*A(k,k-1))/m(k-1);
        g(k) = (b(k)-(g(k-1)*A(k,k-1)))/m(k);
    end

    sol = zeros(1,n);
    sol(n) = g(n);
    for k = n-1:-1:1
        sol(k) = g(k) - (sol(k+1)*A(k,k+1))/m(k);
    end
end
