function [pols] = spline(table)
    
    n = size(table, 1) - 1;

    l = zeros(n+1,1);
    m = zeros(n+1,1);
    d = zeros(n+1,1);
    h = zeros(n,1);

    %for j = 1:n-1
    %    h(j+1) = table(j+1, 1) - table(j, 1);
    %end

    h(1) = table(2, 1) - table(1, 1);
    %h(n) = table(n+1, 1) - table(n, 1);

    for j = 2:n
        h(j) = table(j+1, 1) - table(j, 1);
        l(j) = (h(j))/(h(j-1) + h(j));
        m(j) = 1 - l(j);
        d(j) = (6/(h(j-1) + h(j)))*(((table(j+1, 2) - table(j, 2))/(h(j))) - ...
        ((table(j, 2) - table(j-1, 2))/(h(j-1))));
    end

    A = 2*eye(n+1) + diag(l(1:n), 1)+ diag(m(2:n+1),-1);
    M = tridiag(A, d);
    pols = [];
    for i = 1:n
      aux = [0,0,0,table(i,2)] + ...
      (((table(i+1,2)-table(i,2))/h(i)) - h(i)*((2*M(i)+M(i+1))/6))*[0, 0, 1, -table(i, 1)] + ...
      (M(i)/2)*[0,1,-2*table(i,1), table(i,1)^2] + ...
      ((M(i+1)-M(i))/(6*h(i)))*[1, -3*table(i,1), 3*table(i,1)^2, -table(i,1)^3];
      pols = [pols; aux];
    end
    plotSpline(pols, table);
end

function plotSpline(pols, table)
    n = size(pols, 1);
    for i = 1:n
        r = linspace(table(i,1), table(i+1,1), 100);
        plot(r, polyval(pols(i,:), r));
        hold on;
    end
    hold off;
end

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