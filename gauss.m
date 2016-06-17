function [sol, gaussmat, punt] = gauss(A, b) %Eliminación gaussiana
    n = size(A)(1);
    punt = 1:n;
    
    for i = 1:n-1
        [m pos] = max(abs(A(i:n, i)));
    
        A([i, pos+i-1], :) = A([pos+i-1, i], :);
        punt([i, pos+i-1]) = punt([pos+i-1, i]);
        
        for j = i+1:n
            A(j, i) = (A(j, i)/A(i, i));
        end
        
        for j = i+1:n
            A(j, i+1:n) = A(j, i+1:n) - (A(i, i+1:n)*A(j, i));
        end
    end
    
    sol = gausssolver(A, b, punt);
    gaussmat = A;
end
    
function sol = gausssolver(A, b, punt)
    n = size(A)(1);
    w = zeros(n,1);
    
    w(1) = b(punt(1));
    
    for i = 2:n
        w(i) = b(punt(i)) - dot(A(i, 1:i-1), w(1:i-1));
    end
    
    sol = zeros(n,1);
    
    sol(n) = w(n)/A(n, n);
    
    for i = n-1:-1:1
        sol(i) = (w(i)-dot(A(i, i+1:n), sol(i+1:n)))/A(i, i);
    end
end