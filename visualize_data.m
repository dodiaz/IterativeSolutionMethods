clc
clear all
close all

A = load('A_matrix.txt');

f1 = load('step11957_GS_phi_data.txt');
f2 = load('step21098_PJ_phi_data.txt');
f3 = load('step287_SOR_phi_data.txt');
f4 = load('step34_MG_phi_data.txt');
f5 = load('step113_CG_phi_data.txt');

f1 = f1/64/64;
f2 = f2/64/64;
f3 = f3/64/64;
f4 = f4/64/64;
f5 = f5/64/64;

sum1 = sum(f1,'all');
sum2 = sum(f2,'all');
sum3 = sum(f3,'all');
sum4 = sum(f4,'all');
sum5 = sum(f5,'all');





%% 
error_CG = load('step109_CG_error_data.txt');
error_PJ = load('step21097_PJ_error_data.txt');
error_SOR = load('step304_SOR_error_data.txt');
error_GS = load('step11141_GS_error_data.txt');


figure 
semilogx(error_CG)
xlabel('log(k)')
ylabel('infinity norm of error at step k')
set(gca, 'FontSize', 12)
legend({'Conjugate gradient method'}, 'Location', 'southwest')
figure 
semilogx(error_PJ);
hold on
semilogx(error_GS);
semilogx(error_SOR);
xlabel('log(k)')
ylabel('infinity norm of error at step k')
set(gca, 'FontSize', 12)
legend({'Point Jacobi method', 'Gauss-Seidel', 'Successive over-relaxation (\omega = 1.94)'}, 'Location', 'northeast')
hold off





%% plot the solution

phi = load('step109_CG_phi_data.txt');
phi2 = load('step11141_GS_phi_data.txt');
phi3 = load('step21097_PJ_phi_data.txt');
phi4 = load('step304_SOR_phi_data.txt');



[X,Y] = meshgrid(0:1/63:1);

figure
surf(X,Y,phi)
set(gca, 'FontSize', 11)
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.2 0.4 0.6 0.8 1])
colormap jet
colorbar

