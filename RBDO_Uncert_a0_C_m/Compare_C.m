% compare C plots
clc;clear; close all;

test = dlmread('MC_15000Samples.txt');
C_samples = test(:, 1);
epsF_samples = test(:, 2);

n = 5000 - 1;
[fC, xiC] = ksdensity(C_samples(end-n:end));

fname = 'C:\Users\berkc\Dropbox\Vandy\Research\Optimization\Uniaxial_final\RBDO_April25_Uncert_a0\Figs';
%% Plot
figure

mu_C = 4.5628e-9; cov_C = 0.10;
sigma_C = mu_C * cov_C;

x = linspace (mu_C-6*sigma_C, mu_C+6*sigma_C);
plot(x, normpdf (x,mu_C,sigma_C),'Linewidth',2); hold on;
plot(xiC, fC,'--','Linewidth',2)
xlim([2e-9 7e-9]);

legend('PDF used for optimization', 'PDF from calibration result'...
    ,'location', 'northwest', 'FontName', 'Times New Roman', ...
                    'FontSize',10,'Interpreter', 'LaTeX')
% legend('Crack growth prognosis using Forman model','Test data',...
%     'location', 'northwest')
xlabel('C', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('PDF', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
saveas(gcf, fullfile(fname,'\Ccomparison'), 'pdf') %Save figure 

