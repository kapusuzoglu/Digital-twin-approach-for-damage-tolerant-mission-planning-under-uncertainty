%% Load the test data

% the earlier version starts estimation at 0 cycles, which is not a long...
% crack, this version starts at the first recorded visible cycle.

clear; clc;
close all;

format long
% % read the csv file
M = csvread('032_Marc_27_2019.csv');
n_test = M(1:end, 1);

% % diameter of the initial hole + introduced crack size
% a0 = M(1, 2);

% add diameter
a_test = M(1:end,2);

% % starting point of the estimation
mu_a0 = M(5,2);

mu_C = 4.5628e-9; cov_C = 0.10;

N_Samples = 40; % number of MC samples
U = lhsdesign(N_Samples, 2);

a0_Samples = icdf('normal', U(:, 1), mu_a0, 1e-3); % std = 1e-3 [m]
C_Samples = icdf('normal', U(:, 2), mu_C, mu_C*cov_C);

[a0, C] = ndgrid(a0_Samples,C_Samples);
inp = [a0(:) C(:)];

% 2nd mission
% a_start = 6.096e-3;

% corresponding cycle
n_start = n_test(5, 1);
% n_start = 25000;

% end point of the cycle
nf = n_test(end);

% Array of Cycles ( including 0 cycles ) 
% n_test = M(:, 1);


GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');

GP = {GPStruc_K1 GPStruc_W};

% Fmax = [3634.3 5879.4 4433.9 ...
%         3427.8 5911.2 4429.9 ...
%         3570.5 5827.7 4459.9 ...
%         3600.0 6200.0 4600.0];

% Fmax = [3500.0 5200.0 4500.0 ...
%         3500.0 5200.0 4500.0 ...
%         3500.0 5200.0 4500.0 ...
%         3500.0 5200.0 4500.0];
    
%     Fmax = [3500.3 6100.4 4500.9 ...
%         3500.8 6100.2 4500.9 ...
%         3500.5 6100.7 4500.9 ...
%         3500.0 6100.0 4500.0];
%     Fmax = repmat(5000,1,12);

Fmax = [3500.0 5200.0 4500.0];
Fmin = 0.5*Fmax;

% Nf = n_start + [5213 5213+5247 5213+5247+5325 5213+5247+5325+5310];
% Nf = [1736 1743 1734 1765 1746 1736 1719 1799 1807 1770 1770 1770];
Nf = [1736 1743 1734];



%% Using Forman model to predict crack propagation
tic
% [a_Forman_all, n,Work,W_Sorted_mean,W_Sorted_median] = ...
% Forman_K_gp(GP, Fmax, Fmin, a_start, n_start, n_start+5001, 3.42, ...
%             4.5628e-9, 0.10, N_Samples);
[a_Forman_all, n, W_tot] = ...
Forman_K_gp(GP, Fmax, Fmin, n_start, Nf, 3.42, ...
            inp, N_Samples);
%         4.5628e-9
toc
finalCrackSizes=a_Forman_all(:,end);
a_Forman_mean = mean(a_Forman_all);
finalMeanCrack = a_Forman_mean(end)
a_Forman_std = std(a_Forman_all);
finalStdCrack = a_Forman_std(end)
COV = finalStdCrack/finalMeanCrack;


% fname = 'C:\Users\berkc\Dropbox\Vandy\Research\Optimization\Uniaxial_final\RBDO_April15\Figs';
%% Plot
figure
% ub = max(a_Forman_all);
% lb = min(a_Forman_all);

ub = a_Forman_mean + a_Forman_std;
lb = a_Forman_mean - a_Forman_std;

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
xlabel('Number of Cycles', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Crack Size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\Prognosis_Mission1_RV_a0'), 'pdf') %Save figure               
                
                
% Posterior crack size
[y,x]=ksdensity(finalCrackSizes);
figure; plot(x,y,'k','Linewidth',3)
xlabel('Final crack size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('PDF', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')

set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\PDF_af'), 'pdf') %Save figure


%% Function: Forman_K_gp
% function [a_all,n,Work,W_Sorted_mean,W_Sorted_median] = ...
%          Forman_K_gp(T,Fmax, Fmin, a0, n0, nf, m,...
%          C_mu, C_cov, N_Samples)

function [a_all,n,W_tot] = ...
         Forman_K_gp(T,Fmax, Fmin, n0, nf, m,...
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
% missions = 4;
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

% % 3 blocks for each mission
% W_tot_1 = zeros(N_Samples, missions);
% W_tot_2 = zeros(N_Samples, missions);
% W_tot_3 = zeros(N_Samples, missions);

nf_new = [0 nf];
nf_cum = cumsum(nf_new);
%     ind = 1;

Fmax = repmat(Fmax,N_Samples^2,1);
Fmin = repmat(Fmin,N_Samples^2,1);

C = inp(:,2);
% mu_e = log(1 / sqrt(1 + 0.1^2));
% sigma_e = log(1 + 0.1^2);
tic
% Do for each block loading
for k = 1:length(nf)
    
    % %         for i = (nf - n0)/(3*missions) *(k-1) +1: 1: (nf - n0)/(3*missions) *(k)
    %         % pass the final crack size of the current block
    %         if mod(k,3)==1 && k~=1
    %             a_all(j,ind) = a0(1+ceil(k/4),1);
    %         end
    
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
        dadn = (C .* (deltaK).^m./ ((1 - R(k)) * 67 - deltaK) );%*lognrnd(mu_e,sigma_e);
        
        if (sum(a_dummy>0.025))>=1
            idx = a_dummy>0.025;
            dadn(idx,1)=0;
        end
        
%         end
        
        a_all(:, nf_cum(k) + i + 1 ) = a_dummy + 2 * dadn;
        
%         W2_k(j,i,k) = predict(GPR_W, [Fmax(k) 1000*a_all(j, i+1)]);
%         % work done from 0 to Fmin
%         W1_k(j,i,k) = predict(GPR_W, [Fmin(k) 1000*a_all(j, i+1)]);
%         % Work done = W2-W1
%         W_all(j,i,k) = 0.5*(W2_k(j,i,k)-W1_k(j,i,k));
        
        % Predict work done in the next cycle (only for the increase in
        % force)
        % work done from 0 to Fmax
        % convert cracksize from micrometers to meters
        W2_k(:,i,k) = predict(GPR_W, [Fmax(:,k) a_all(:, i+1)]);
        % work done from 0 to Fmin
        W1_k(:,i,k) = predict(GPR_W, [Fmin(:,k) a_all(:, i+1)]);
        % Work done = W2-W1
        W_all(:,i,k) = 0.5*(W2_k(:,i,k)-W1_k(:,i,k));
    end
    %         ind = nf_cum(k) + i + 1;
    a_all(:, nf_cum(k) + i + 1 ) = a_all(:, nf_cum(k) + i + 1 );% + normrnd(0,8.678052985519797e-04);
    if mod(k,3) == 0
        fprintf('Mission %d is completed.\n', round(k/3))
        toc
    end
end

W_tot_1 = sum(W_all(:,:,1),2);
W_tot_2 = sum(W_all(:,:,2),2);
W_tot_3 = sum(W_all(:,:,3),2);

% pass mean of work dones of all samples
W_tot = [mean(W_tot_1); mean(W_tot_2); mean(W_tot_3)];

end