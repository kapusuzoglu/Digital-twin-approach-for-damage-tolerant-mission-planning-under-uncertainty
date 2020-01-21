% plot KL vs. # of particles

% in mm
KL_10avg = (0.006107167061003e3 + 0.006159492737520e3 + 0.005924915851059e3 + 0.006168050510559e3)/4;
KL_15avg = 0.006315534016836e3;
KL_20avg = 0.006135681816389e3;
KL_25avg = 0.006155449507640e3;
KL_30avg = 0.006149442408928e3;
KL_35avg = 0.006147396277373e3;
KL_40avg = 0.006127821907048e3;
KL_50avg = 0.006141384311096e3;
KL_60avg = 0.006141055417478e3;
KL_70avg = 0.006135253885457e3;
KL_80avg = 0.006167653493403e3;
KL_90avg = 0.006166605360706e3;
KL_100avg = 0.006147688946223e3;

fname = 'C:\Users\berkc\Dropbox\Vandy\Research\Optimization\Uniaxial_final\RBDO_April15\Figs';
%% Figures
figure
x = [10^2 15^2 20^2 25^2 30^2 35^2 40^2 50^2 60^2 70^2 80^2 90^2 100^2];
y = [KL_10avg  KL_15avg KL_20avg KL_25avg KL_30avg KL_35avg KL_40avg KL_50avg KL_60avg KL_70avg KL_80avg KL_90avg KL_100avg];
plot(x,y,'-ks','LineWidth',2);
xlabel('Number of samples', 'FontName', 'Times New Roman', ...
                    'FontSize',14,'Color','k', 'Interpreter', 'LaTeX');
ylabel('Mean of final crack size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',14,'Color','k', 'Interpreter', 'LaTeX');
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\Prognosis_Mission1_RV_a0_MEAN'), 'pdf') %Save figure 

% in mm                
KL_10avg = (0.001485160925523e3 + 0.001649266837056e3 + 0.001495649540744e3 + 0.001574656484317e3)/4;
KL_15avg = 0.001665599243518e3;
KL_20avg = 0.001343882465303e3;
KL_25avg = 0.001352074965928e3;
KL_30avg = 0.001363503305993e3;
KL_35avg = 0.001422532499497e3;
KL_40avg = 0.001450893752040e3;
KL_50avg = 0.001373506551629e3;
KL_60avg = 0.001386911172401e3;
KL_70avg = 0.001422745876643e3;
KL_80avg = 0.001421691763795e3;
KL_90avg = 0.001416580689387e3;
KL_100avg = 0.001408221873577e3;

%% Figures
figure
x = [10^2 15^2 20^2 25^2 30^2 35^2 40^2 50^2 60^2 70^2 80^2 90^2 100^2];
y = [KL_10avg KL_15avg KL_20avg KL_25avg KL_30avg KL_35avg KL_40avg KL_50avg KL_60avg KL_70avg KL_80avg KL_90avg KL_100avg];
plot(x,y,'-ks','LineWidth',2);
xlabel('Number of samples', 'FontName', 'Times New Roman', ...
                    'FontSize',14,'Color','k', 'Interpreter', 'LaTeX');
ylabel('Standard deviation of final crack size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',14,'Color','k', 'Interpreter', 'LaTeX');
                
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\Prognosis_Mission1_RV_a0_STD'), 'pdf') %Save figure 