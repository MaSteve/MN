function funs = cholesky()
    funs.choleskyM=@choleskyM;
    funs.choleskysolver=@choleskysolver;
end

function [sol, cholmat] = choleskyM(A, b) %Resuelve sistemas por el método de Cholesky.
    n = size(A)(1);
    ok = true;

    for i = 1:n
        A(i,i) = (A(i,i) - dot(A(i, 1:i-1), A(i, 1:i-1)));
        if A(i,i) < 0
            disp('No se puede hacer Cholesky')
            ok = false;
            break;
        end
        A(i,i) = sqrt(A(i,i));

        for j = i+1:n
            A(j,i) = (A(i,j) - dot(A(i, 1:i-1),A(j, 1:i-1)))/A(i,i);
        end
    end
    if ok
        sol = choleskysolver(A, b);
        cholmat = A;
    end
end

function sol = choleskysolver(Amod, b) %Permite resolver varios sistemas con la misma matriz
    n = size(Amod)(1);

    w = zeros(n,1);
    for i = 1:n
        w(i) = (b(i) - dot(Amod(i, 1:i-1), w(1:i-1)))/Amod(i, i);
    end

    sol = zeros(n,1);
    for i = n:-1:1
        sol(i) = (w(i) - dot(Amod(i+1:n, i), sol(i+1:n)))/Amod(i, i);
    end
end
