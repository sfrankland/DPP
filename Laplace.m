function d = Laplace(A)
[~,n] = size(A);
if n == 1
    d = A;
else
    d = 0;
    for j = 1:n
        A1 = A(2:n, [1:j-1, j+1:n]);
        d = d+(-1)^(j+1)*A(1,j)*Laplace(A1);
    end
end
