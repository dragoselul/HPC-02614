% matrix times vector
%
n = 500:100:5000;

j = 1;
for i = n
    A = ones(i);
    b = ones(i,1);
    c = @() A*b;
    tmv(j) = timeit(c);
    fprintf("%5d: %f\n", i, tmv(j));
    j = j + 1;
end

pmv = (n.^2) ./ tmv * 2.0e-6;

figure
subplot(2,1,1);
plot(n,tmv, 'x-')
title('matrix times vector')
ylabel('runtime')
subplot(2,1,2);
plot(n,pmv, 'o-')
ylabel('performance')
xlabel('n')

