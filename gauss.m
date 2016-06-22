function funs = gauss()
    funs.gaussM=@gaussM;
    funs.gausssolver=@gausssolver;
end

function [sol, gaussmat, punt] = gaussM(A, b) %Eliminación gaussiana
    n = size(A)(1);
    punt = 1:n;

    for i = 1:n-1
        [m pos] = max(abs(A(punt(i:n), i)));

        %A([i, pos+i-1], :) = A([pos+i-1, i], :); %no permutar
        punt([i, pos+i-1]) = punt([pos+i-1, i]);

        for j = i+1:n
            A(punt(j), i) = (A(punt(j), i)/A(punt(i), i));
        end

        for j = i+1:n
            A(punt(j), i+1:n) = A(punt(j), i+1:n) - (A(punt(i), i+1:n)*A(punt(j), i));
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
        w(i) = b(punt(i)) - dot(A(punt(i), 1:i-1), w(1:i-1));
    end

    sol = zeros(n,1);

    sol(n) = w(n)/A(punt(n), n);

    for i = n-1:-1:1
        sol(i) = (w(i)-dot(A(punt(i), i+1:n), sol(i+1:n)))/A(punt(i), i);
    end
end
