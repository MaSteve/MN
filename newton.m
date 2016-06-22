function funs = newton()
    funs.newtonP=@newtonP;
    funs.addPoint=@addPoint;
    funs.fromFunction=@fromFunction;
end

function [pol, Pi, tab] = newtonP(table) %Vector de dos columnas
    n = size(table)(1) - 1;

    pol = [];
    Pi = [1];

    for k = 0:n-1
        pol = [0, pol] + Pi*table(1,2);
        Pi = [Pi, 0] - [0, Pi.*table(k+1,1)];
        for i = 1:n-k
            table(i,2) = (table(i,2) - table(i+1,2))/(table(i,1) - table(i+k+1,1));
        end
    end

    pol = [0, pol] + Pi*table(1,2);
    Pi = [Pi, 0] - [0, Pi.*table(n+1,1)];
    tab = table;

    plotNewton(pol, table);
end

function plotNewton(pol, table)
    r = linspace(min(table(:, 1)), max(table(:, 1)), 100);
    plot(r, polyval(pol, r));
end

function [pol, Pi, tab] = fromFunction(f, v) %vector de puntos
	n = size(v, 2);
    f = eval(['@(x)' f]);
    table = zeros(n, 2);
    for i = 1:n
        table(i,1) = v(i);
        table(i,2) = f(table(i,1));
    end
    [pol, Pi, tab] = newtonP(table);
    hold on;
    r = linspace(table(1,1), table(n,1), 100);
    plot(r, f(r(:)));
    hold off;
end

function [npol, nPi, tab] = addPoint(val, img, pol, Pi, table)
    n = size(table)(1);
    table = [table; [val, img]];

    for i = n:-1:1
        table(i,2) = (table(i,2) - table(i+1,2))/(table(i,1) - table(n+1,1));
    end

    pol = [0, pol] + Pi*table(1,2);
    Pi = [Pi, 0] - [0, Pi.*table(n+1,1)];
    npol = pol;
    nPi = Pi;
    tab = table;

    plotNewton(pol, table);
end
