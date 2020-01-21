%% Load the test data

% the earlier version starts estimation at 0 cycles, which is not a long...
% crack, this version starts at the first recorded visible cycle.

clear; clc;
close all;

% % read the csv file
% M = csvread('032_Apr_25_2019.csv');
% n_test = M(1:end, 1);

% % diameter of the initial hole + introduced crack size
% a0 = M(1, 2);
% 
% add diameter
% a_test = M(1:end,2);

N_Samples = 10; % number of MC samples
U = lhsdesign(N_Samples, 6);

% Diagnosis w/ uncertainty
% % hybrid with a0 + C
% a_start = [4.572e-3; 5.59e-3; 9.906e-3; 13.462e-3; 16.256e-3];
% only RDO with a0 + C
a_start = [4.826e-3; 6.096e-3; 10.668e-3; 14.986e-3; 18.288e-3];
a0_Diagnosis1 = norminv(U(:, 1), a_start(1,1), 1e-3); % std = 1e-3 [m]
a0_Diagnosis2 = norminv(U(:, 2), a_start(2,1), 1e-3); % std = 1e-3 [m]
a0_Diagnosis3 = norminv(U(:, 3), a_start(3,1), 1e-3); % std = 1e-3 [m]
a0_Diagnosis4 = norminv(U(:, 4), a_start(4,1), 1e-3); % std = 1e-3 [m]

a0_Diagnosis = [a0_Diagnosis1 a0_Diagnosis2 a0_Diagnosis3 a0_Diagnosis4];

% % starting point of the estimation
mu_a0 = a_start(1,1);

mu_C = 4.5628e-9; cov_C = 0.10;

a0_Samples = icdf('normal', U(:, 5), mu_a0, 1e-3); % std = 1e-3 [m]
C_Samples = icdf('normal', U(:, 6), mu_C, mu_C*cov_C);

[a0, C] = ndgrid(a0_Samples,C_Samples);
inp = [a0(:) C(:)];

% Diagnosis
[a0_diag1, C] = ndgrid(a0_Diagnosis(:,1),C_Samples);
[a0_diag2, C] = ndgrid(a0_Diagnosis(:,2),C_Samples);
[a0_diag3, C] = ndgrid(a0_Diagnosis(:,3),C_Samples);
[a0_diag4, C] = ndgrid(a0_Diagnosis(:,4),C_Samples);
a0_diag = [a0_diag1(:) a0_diag2(:) a0_diag3(:) a0_diag4(:)];

% corresponding cycle
n_start = 25000 + 1250;

% end point of the cycle
% nf = n_test(end);

% Array of Cycles ( including 0 cycles ) 
% n_test = M(:, 1);

GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');

GP = {GPStruc_K1,GPStruc_W};

% % RBDO + RDO only C
% Fmax = [3411.1 5000.2 4000.0 3445.3 5000.0 4008.7 3690.6 5000.0 4000.0 3500 5200 4500];
% Nf = [1179 4034 1125 792 2676 751 1528 4912 1528 422 1266 422];

% % RBDO + RDO for C + a0
% Fmax = [3263.2 5000.2 4089.5 3191.4 5000.0 4093.8 3000 5000.0 4188.9 3500 5200 4500];
% Nf = [1238 4115 1203 819 2713 759 1545 5260 1654 422 1266 422];

% % RDO only C
% Fmax = [3411.1 5000.2 4000.0 3445.3 5000.0 4008.7 3000.0 5000.0 4000.0 3022 5200 4005.9];
% Nf = [1179 4034 1125 792 2676 751 2015 4438 1999 417 1306 389];

% hybrid for C + a0
Fmax = [3000 5000 4004.3 3000 5000.0 4000 3070.9 5000.0 4000.2 3000.7 5000 4000];
Nf = [1322 1745 1388 2659 3059 2458 1734 1698 1800 1044 1312 895];

% only RDO for C + a0
% Fmax = [3000.1 5027.1 4076.9 3131.2 5000.0 4000 3063.6 5000 4000 3000 5000 4000];
% Nf = [1338 1766 1255 2610 3272 2383 1741 1703 1800 1038 1301 962];
n_test = n_start + [0, sum(Nf(1,1:3)), sum(Nf(1,1:6)), sum(Nf(1,1:9)), sum(Nf(1,1:12))];

Fmin = 0.5*Fmax;

mu_C = 4.5628e-9; cov_C = 0.10;
U = lhsdesign(N_Samples, 1);
C_Samples = icdf('normal', U(:, 1), mu_C, mu_C*cov_C);


%% Using Forman model to predict crack propagation
tic
% [a_Forman_all, n,Work,W_Sorted_mean,W_Sorted_median] = ...
% Forman_K_gp(GP, Fmax, Fmin, a_start, n_start, n_start+5001, 3.42, ...
%             4.5628e-9, 0.10, N_Samples);
[a_Forman_all, n,Work,W_Sorted_mean,W_Sorted_median] = ...
Forman_K_gp(GP, Fmax, Fmin, a0_diag,n_start, Nf, 3.42, ...
            inp, N_Samples);
toc
finalCrackSizes=a_Forman_all(:,end);
a_Forman_mean = mean(a_Forman_all);
finalMeanCrack = a_Forman_mean(end)
a_Forman_std = std(a_Forman_all);


a_test = [a_start'];

fname = 'C:\Users\berkc\Dropbox\Vandy\Research\Optimization\Uniaxial_final\RBDO_April25_Uncert_a0\Figs';
%% Plot
figure

% ub = max(a_Forman_all);
% lb = min(a_Forman_all);

ub = a_Forman_mean + 2*a_Forman_std;
lb = a_Forman_mean - 2*a_Forman_std;

xconf = [n n(end:-1:1)] ;         
yconf = [ub*1e3 lb(end:-1:1)*1e3];

p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.5];      
p.EdgeColor = 'none'; 
hold on;

% plot(n, a_Forman_all, '-k', 'linewidth', 2);hold on;
plot(n, a_Forman_mean*1e3, '-k', 'linewidth', 2);hold on;
plot(n_test, a_test*1e3, '-bo')

legend('Uncertainty Bounds ($\mu \pm 2\sigma$)', 'Model prediction'...
    ,'Test data','location', 'northwest', 'FontName', 'Times New Roman', ...
                    'FontSize',10,'Interpreter', 'LaTeX')
% legend('Crack growth prognosis using Forman model','Test data',...
%     'location', 'northwest')
xlabel('Number of Cycles', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Crack Size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\onlyRDO_Prognosis_Optimization_C_and_a0_withDiagnosis'), 'pdf') %Save figure 


%% Function: Forman_K_gp
% function [a_all,n,Work,W_Sorted_mean,W_Sorted_median] = ...
%          Forman_K_gp(T,Fmax, Fmin, a0, n0, nf, m,...
%          C_mu, C_cov, N_Samples)

function [a_all,n,Work,W_Sorted_mean,W_Sorted_median] = ...
         Forman_K_gp(T,Fmax, Fmin, a0,n0, nf, m,...
         inp, N_Samples)
     
% ------ fitrgp ------
% load('gprMdl_K1','gprMdl_K1')
GPR_K1 = T{1}.('gprMdl_K1');
GPR_W = T{2}.('gprMdl_W');
% ------ fitrgp ------
                         
% ------ Physical Parameters ------
% Strss Ratio
R = Fmin./Fmax;

% 3 missions for 3 manevours
missions = 4;
maneuver = 3;
% Fmax = repmat(Fmax,1,missions);
% Fmin = repmat(Fmin,1,missions);
% R = repmat(R,1,missions);

% nf = n0+(nf-n0);

% cycles counts at which the crack size needs to be estimated
n = n0:1:n0+sum(nf);


% place holders for estimated crack size
% each row is a trace of crack growth corresponding to a set of MC Samples
a_all = zeros(N_Samples^2, length(n));

% delta_K = zeros(N_Samples, length(n)-1);
% da_dn = zeros(N_Samples, length(n)-1);
a_all(:, 1) = inp(:,1);

% Work done
W1_k = zeros(N_Samples^2, length(n)-1,length(Fmax));
W2_k = zeros(N_Samples^2, length(n)-1,length(Fmax));
W_all = zeros(N_Samples^2, length(n)-1,length(Fmax));


% 3 blocks for each mission
W_tot_1 = zeros(N_Samples^2, missions);
W_tot_2 = zeros(N_Samples^2, missions);
W_tot_3 = zeros(N_Samples^2, missions);

nf_new = [0 nf];
nf_cum = cumsum(nf_new);
%     ind = 1;

Fmax = repmat(Fmax,N_Samples^2,1);
Fmin = repmat(Fmin,N_Samples^2,1);
format long
C = inp(:,2);

ind = 1;  
% Do for each block loading
for k = 1:length(nf)
  
    % pass the final crack size of the current block
    if mod(k,3)==1 && k~=1
        mean(a_all(:,ind))
        a_all(:,ind) = a0(:,1+ceil(k/4),1);
    end
    
    i = 0;
    
    while i < nf_new(k+1)
        i=i+1;
        
        a_dummy = a_all(:, nf_cum(k) + i);
        
        % gp model input requirement: force in lbs, crack length in mm
        [Kmax_mean] = predict(GPR_K1, [Fmax(:,k) 1000*a_dummy]);
        [Kmin_mean] = predict(GPR_K1, [Fmin(:,k) 1000*a_dummy]);
        
        deltaK = Kmax_mean - Kmin_mean;
        % deltaK from gp model use MPa*sqrt(m) in the Forman Law
        
        deltaK = deltaK/(1e6);
        
        % determine da/dn value - a here is HALF crack length
        %                 dadn = (C * (deltaK)^m/ ((1 - R(k)) * 67 - deltaK) );% + normrnd(6.505792002037899e-08,1.986698609510128e-08); %+ eps_Forman);% + abs(eps_obs);
        dadn = (C .* (deltaK).^m./ ((1 - R(k)) * 67 - deltaK) );
            
        if (sum(a_dummy>0.035))>=1
            idx = a_dummy>0.035;
            dadn(idx,1)=0;
        end
        
        a_all(:, nf_cum(k) + i + 1 ) = a_dummy + 2 * dadn;
        
                    W2_k(:,i,k) = predict(GPR_W, [Fmax(:,k) 1000*a_all(:, nf_cum(k) + i + 1 )]);
                    % work done from 0 to Fmin
                    W1_k(:,i,k) = predict(GPR_W, [Fmin(:,k) 1000*a_all(:, nf_cum(k) + i + 1 )]);
                    % Work done = W2-W1
                    W_all(:,i,k) = 0.5*(W2_k(:,i,k)-W1_k(:,i,k));
    end
    ind = nf_cum(k) + i + 1;
%     a_all(:, nf_cum(k) + i + 1 ) = a_all(:, nf_cum(k) + i + 1 );% + normrnd(0,8.678052985519797e-04);
    
end

% Total amount of work done for each blocks and for different samples
% for each mission
kk = 0; jj=0;
for ii=1:missions
    jj = jj+1;
    W_tot_1(:,jj) = sum(W_all(:,:,1+kk),2);
    W_tot_2(:,jj) = sum(W_all(:,:,2+kk),2);
    W_tot_3(:,jj) = sum(W_all(:,:,3+kk),2);
    kk = kk + maneuver;
end
Work = [W_tot_1 W_tot_2 W_tot_3];
% pass mean of work dones of all samples
W_tot_mean = mean(Work);

W_tot_median = median(Work);

% W_tott=Work
% W_tot = [W_tott(1,1), W_tott(1,5), W_tott(1,9),W_tott(1,2), W_tott(1,6),...
%     W_tott(1,10), W_tott(1,3), W_tott(1,7), W_tott(1,11),...
%     W_tott(1,4), W_tott(1,8), W_tott(1,12)];

% W_Sorted = [W_tott(1,1), W_tott(1,4), W_tott(1,7),W_tott(1,2), W_tott(1,5),...
%     W_tott(1,8), W_tott(1,3), W_tott(1,6), W_tott(1,9)];
W_Sorted_mean = [W_tot_mean(1,1), W_tot_mean(1,2), W_tot_mean(1,3)];
W_Sorted_median = [W_tot_median(1,1), W_tot_median(1,2), W_tot_median(1,3)];

end