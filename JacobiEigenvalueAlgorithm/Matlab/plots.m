x = [ 6 16 64 128 256 512 712];
qr = [0.079777 	0.037322	0.52578 	1.97342 	7.62218 	28.4504 	55.895];
dc = [0.003732 	0.011197	0.125963	0.359695	1.39166 	5.19435 	9.41272];
%%
% e-6

y12 = [0.568235  	3.30257 	135.075	916.581 	8369.91 	76915.5 	214435]; 
y22 = [0.343833 	3.68326 	166.88 	1247.22	10393	77948.4 	204061];
y32 = [0.363894 	2.16983 	88.4752 	599.748 	5667.41 	46571.3 	136608];
y42 = [2.47355 	2.87291	79.9158 	456.629 	4680.92 	32113	78365.3];
figure(1)
plot(x, y12, '-o');
hold on 
plot(x,y22, '-o')
hold on 
plot(x,y32, '-o')
hold on 
plot(x,y42, '-o')
hold on 
plot(x,qr, '-*')
hold on 
plot(x,dc, '-*')
grid on
xlabel('n');
ylabel('time elapsed (ms)');
legend('actual', 'Ideal');
title('Time elapsed: e-6')
legend('jacobi-p1', 'jacobi-p2', 'jacobi-p4', 'jacobi-p8', 'qr', 'dc')
hold off
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
%%
% e-3
y1 = [0.186612	2.07699 	74.2969	704.73 	5592.5 	46814.8	125285]; 
y2 = [0.312109 	1.83953	49.1695 	467.021 	4007.66 	28766.1 	78396.3];
y3 = [0.439938 	1.45324	40.1454 	355.453	3242.93 	23425.3	60743.7 ];
y4 = [1.257	1.74202	49.45181	377.366	2616.35 	17287.6 	47560.3];

figure(2)
plot(x, y1, '-o')
hold on 
plot(x,y2, '-o')
hold on 
plot(x,y3, '-o')
hold on 
plot(x,y4, '-o')
hold on 
plot(x,qr, '-*')
hold on 
plot(x,dc, '-*')
grid on
xlabel('n');
ylabel('time elapsed (ms)');
title('Time elapsed: e-3')
legend('jacobi-p1', 'jacobi-p2', 'jacobi-p4', 'jacobi-p8', 'qr', 'dc')
hold off
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));

speed_up2 = y1./y2;
speed_up3 = y1./y3;
speed_up4 = y1./y4;

eff2 = speed_up2./2;
eff3 = speed_up2./4;
eff4 = speed_up2./8;

figure(3)
plot(x, speed_up2, '-o')
hold on 
plot(x, speed_up3, '-o')
hold on 
plot(x, speed_up4, '-o')
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
plot(x, eff3, '-o')
hold on 
plot(x, eff4, '-o')
hold on 
grid on
xlabel('n');
ylabel('Efficiency');
title('Efficiency: e-6')
legend('efficiency-2p', 'efficiency-4p', 'efficiency-8p')
hold off

work712(1)=y1(7);
work712(2)=y2(7);
work712(3)=y3(7);
work712(4)=y4(7);
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