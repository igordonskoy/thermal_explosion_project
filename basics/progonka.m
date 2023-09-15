function X = progonka(A,B)

% function X = progonka(A,B)
% Решение трехдиагональной СЛАУ методом прогонки
if size(A,1) ~= size(B,1)
    error(['Dimensions of matrix: ', num2str(size(A)), '; dimensions of vector: ', num2str(size(B))])
end


for j = 1:size(B,2)
    Bj = B(:,j);
    N = length(Bj); a = zeros(N,1); b = zeros(N,1); a(2) = -A(1,2)/A(1,1); b(2) = Bj(1)/A(1,1);
    for i = 2:(N-1)
        p =(A(i,i-1)*a(i)+A(i,i)); a(i+1) = -A(i,i+1)/p; b(i+1) = (Bj(i) - A(i,i-1)*b(i))/p;
        if p == 0
            error('The matrix is uninvertible');
        end
    end
    Xj = zeros(N,1); Xj(N) = (Bj(N) - A(N,N-1)*b(N))/(A(N,N) + A(N,N-1)*a(N));
    for i = (N-1):(-1):1
        Xj(i) = a(i+1)*Xj(i+1) + b(i+1);
    end
    X(:,j) = Xj;
end