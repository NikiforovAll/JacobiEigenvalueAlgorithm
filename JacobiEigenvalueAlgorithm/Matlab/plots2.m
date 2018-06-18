x = [16 64 128 256 512 712];


bm1 = csvread('.\test-data\bisection_test-e6-th1.csv');
bm2 = csvread('.\test-data\bisection_test-e6-th2.csv');
bm4 = csvread('.\test-data\bisection_test-e6-th4.csv');
bm8 = csvread('.\test-data\bisection_test-e6-th8.csv');
bm_lp = csvread('.\test-data\sstebz_lapacktest.csv');
b1 = bm1(:,2);
b2 = bm2(:,2);
b4 = bm4(:,2);
b8 = bm8(:,2);
bl = bm_lp(:,2);
figure(1)
plot(x, b1, '-o');
hold on 
plot(x, b2, '-o')
hold on 
plot(x,b4, '-o')
hold on 
plot(x,b8, '-o')
hold on 
plot(x,bl, '-*')
grid on
xlabel('n');
ylabel('time elapsed (ms)');
title('Time elapsed: e-6')
legend('bisection-p1', 'bisection-p2', 'bisection-p4', 'bisection-p8', 'sstebz-lapack')
hold off
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));

speed_up2 = b1./b2;
speed_up4 = b1./b4;
speed_up8 = b1./b8;

eff2 = speed_up2./2;
eff4 = speed_up4./4;
eff8 = speed_up8./8;

figure(3)
plot(x, speed_up2, '-o')
hold on 
plot(x, speed_up4, '-o')
hold on 
plot(x, speed_up8, '-o')
hold on 
grid on
xlabel('n');
ylabel('Speedup');
title('Speedup: e-6')
legend('speedup-2p', 'speedup-4p', 'speedup-8p')
hold off

figure(4)
plot(x, eff2, '-o')
hold on 
plot(x, eff4, '-o')
hold on 
plot(x, eff8, '-o')
hold on 
grid on
xlabel('n');
ylabel('Efficiency');
title('Efficiency: e-6')
legend('efficiency-2p', 'efficiency-4p', 'efficiency-8p')
hold off

work712(1)=b1(6);
work712(2)=b2(6);
work712(3)=b4(6);
work712(4)=b8(6);
x_workers = [1 2 4 8];
speed_up_workers = work712(1)./work712;

figure(5)
A = [1, 2];
B = [speed_up_workers(1), speed_up_workers(2)];
plot(x_workers , speed_up_workers, '-o', 'Color', 'b','MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5], 'MarkerSize', 7);

xlim = get(gca,'XLim');
m = (B(2)-B(1))/(A(2)-A(1));
n = B(2) - A(2)*m;
yp1 = m*xlim(1) + n;
yp2 = m*xlim(2) + n;
hold on
grid on
line([xlim(1) xlim(2)],[yp1 yp2], 'Color', 'r')
xlabel('number of threads');
ylabel('Speedup ration');
legend('actual', 'Ideal');
title('Speedup ration: n = 712; eps = e-6;')
figure(6)

speed_up_workers_eff = speed_up_workers./x_workers;

plot(x_workers, speed_up_workers_eff,'-o', 'Color', 'b','MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5], 'MarkerSize', 7);
ylim([0, 1.1]);
hold on 
grid on
line([1 8],[1 1], 'Color', 'r')
xlabel('number of threads');
ylabel('Parallel Efficiency');
legend('actual', 'Ideal');
title('Parallel Efficiency: n = 712; eps = e-6;')