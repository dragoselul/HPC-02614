% matrix times matrix
%
n = 500:500:5000;

j = 1;
for i = n
    A = ones(i);
    tic
    B = A*A;
    tmm(j) = toc;
    fprintf("%5d: %f\n", i, tmm(j));
    j = j + 1;
end

pmm = (n.^3) ./ tmm * 2.0e-6;

figure
subplot(2,1,1);
plot(n,tmm, 'x-')
title('matrix times matrix')
ylabel('runtime')
subplot(2,1,2);
plot(n,pmm, 'o-')
ylabel('performance')
xlabel('n')

