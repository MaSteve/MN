function [sol, lumat] = lufact(A, b) %Resuelve sistemas por el método de factorización LU.
    n = size(A)(1);
    ok = true;
    %L = eye(n)
    %U = zeros(n,n)
    for i = 1:n
        A(i,i) = A(i,i) - (A(i, 1:i-1)*A(1:i-1, i));
        if A(i,i) == 0
            disp('No se puede hacer LU')
            ok = false;
            break;
        end
        %U(i,i) = A(i,i);
        for j = i+1:n
            A(i,j) = A(i,j) - (A(i, 1:i-1)*A(1:i-1, j));
            %U(i,j) = A(i,j);
        end
        for j = i+1:n
            A(j,i) = (A(j,i) - (A(j, 1:i-1)*A(1:i-1, i)))/A(i,i);
            %L(j,i) = A(j,i);
        end
    end
    if ok
        sol = lusolver(A, b);
        lumat = A;
    end
end 
    
function sol = lusolver(Amod, b) %Permite resolver varios sistemas con la misma matriz
    n = size(Amod)(1);
    %Solución sistema triangular inferior con 1's en la diagonal
    w = zeros(n,1);
    for i = 1:n
        w(i) = b(i) - dot(Amod(i, 1:i-1), w(1:i-1));
    end
    
    %L*w
    
    %Solución sistema triangular superior arbitraria
    sol = zeros(n,1);
    for i = n:-1:1
        sol(i) = (w(i) - dot(Amod(i, i+1:n), sol(i+1:n)))/Amod(i, i);
    end
end