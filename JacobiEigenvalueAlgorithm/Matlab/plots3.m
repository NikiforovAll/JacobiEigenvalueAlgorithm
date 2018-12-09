x = [16 64 128 256 512 712];
shift = 2;
x = x(1:1,1:shift);
prefix = '.\test-data2\04-12-18-output-';
file1 = strcat(prefix, 'bisection_test.csv');
file2 = strcat(prefix, 'sstebz_lapacktest.csv');
file3 = strcat(prefix, 'ssteqr_lapacktest.csv');
file4 = strcat(prefix, 'stedc_lapacktest.csv');
file5 = strcat(prefix, 'parallel_jacob_musictest.csv');
% horrible code
bm1 = csvread(file1);
b1 = bm1(1:shift,2);
bm2 = csvread(file2);
b2 = bm2(1:shift,2);
bm3 = csvread(file3);
b3 = bm3(1:shift,2);
bm4 = csvread(file4);
b4 = bm4(1:shift,2);
bm5 = csvread(file5);
b5 = bm5(1:shift,2);

figure(1)
plot(x, b1, '-o');
hold on 
plot(x, b2, '-o')
hold on 
plot(x, b3, '-o')
hold on 
plot(x, b4, '-o')
hold on 
plot(x, b5, '-o')

grid on
xlabel('n');
curtick = get(gca, 'YTick');
ylabel('time elapsed (ms)');
title('Time elapsed: e-6')
hold off
legend('bisection-th8', 'sstebz', 'ssteqr', 'stedc', 'parallel-jacobi')
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
