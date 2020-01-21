%% Load the test data

% Estimate the initial crack sizes of 2nd, 3rd and 4th missions
% by using the results from the optimization, i.e. loading,cycles, and
% final mean crack sizes of previous missions

clear; clc;
close all;
format long
% % read the csv file
M = csvread('032_Marc_27_2019.csv');
n_test = M(1:end, 1);

% diameter of the initial hole + introduced crack size
% a0 = M(1, 2);

% add diameter
a_test = M(1:end,2);

% % starting point of the estimation
% a_start = a_test(5:end,1); % the initial crack size of 1st mission
% a_start = 6.666292293350e-3;   % the initial crack size of 2nd mission
a_start = 8.590680307865e-3; % the initial crack size of 3rd mission

% mu_a0 = a_start; cov_a0 = 0.1;
% a_start = icdf('Normal', lhsdesign(1,1), mu_a0, mu_a0*cov_a0);

% corresponding cycle
n_start = n_test(5, 1); % 1st mission
% n_start = n_test(6, 1); % 2nd mission
% n_start = n_test(7, 1); % 3rd mission

% end point of the cycle
nf = n_test(end);

% Array of Cycles ( including 0 cycles ) 
% n_test = M(:, 1);


GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');

GP = {GPStruc_K1 GPStruc_W};

% Fmax = [3411.1 5000.2 4000.0]; % 1st mission
% Fmax = [3445.3 5000.0 4008.7]; % 2nd mission
Fmax = [3000.0 5000.0 4000.0]; % 3rd mission
Fmin = 0.5*Fmax;


% Nf = [1179 4034 1125]; % 1st mission
% Nf = [792 2676 751]; % 2nd mission
Nf = [1350 4050 1350]; % 3rd mission

n_start = n_start + 792 + 2676 + 751

N_Samples = 2; % number of MC samples


mu_C = 4.7e-9; cov_C = 0.001;
U = lhsdesign(N_Samples, 1);
C_Samples = icdf('normal', U(:, 1), mu_C, mu_C*cov_C)


%% Using Forman model to predict crack propagation
tic
% [a_Forman_all, n,Work,W_Sorted_mean,W_Sorted_median] = ...
% Forman_K_gp(GP, Fmax, Fmin, a_start, n_start, n_start+5001, 3.42, ...
%             4.5628e-9, 0.10, N_Samples);
[a_Forman_all, n, W_tot] = ...
Forman_K_gp(GP, Fmax, Fmin, a_start, n_start, Nf, 3.42, ...
            C_Samples, N_Samples);
%         4.5628e-9
toc
finalCrackSizes=a_Forman_all(:,end)
a_Forman_mean = mean(a_Forman_all);
finalMeanCrack = a_Forman_mean(end)
a_Forman_std = std(a_Forman_all);

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
% legend('Crack growth prognosis using Forman model','Test data',...
%     'location', 'northwest')
xlabel('Number of Cycles', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Crack Size [\textit{mm}]', 'FontName', 'Times New Roman', ...
                    'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')

                
% a_F = finalCrackSizes * 1e3; % in [mm]
% a_crit = a_crit*1e3; % in [mm]

% figure();
% [fm,xm] = ksdensity(a_F(:));
% [~, index] = min(abs(xm-a_crit));
% plot(xm,fm,'.-r')
% hold on
% line([a_crit a_crit], [min(fm) max(fm)*(1+0.1) ]);
% text(a_crit+a_crit*.00003,max(fm),'\rightarrow a_{crit}','FontSize',16);
% hold on
% p = xm(index:end);
% q = fm(index:end);
% H = area(p,q);
% set(H(1),'FaceColor',[1 0.5 0])
% xlabel('Crack size [\textit{mm}]', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('PDF', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% 
%                 
% figure();
% [fm,xm] = ksdensity(Work(:,1));
% [~, index] = min(abs(xm-W_Sorted_median(1,1)));
% plot(xm,fm,'.-r')
% hold on
% line([W_Sorted_median(1,1) W_Sorted_median(1,1)], [min(fm) max(fm)*(1+0.1) ]);
% % text(W_Sorted_median(1,1)+W_Sorted_median(1,1)*.00003,max(fm),'\rightarrow a_{crit}','FontSize',16);
% hold on
% p = xm(index:end);
% q = fm(index:end);
% H = area(p,q);
% set(H(1),'FaceColor',[1 0.5 0])
% xlabel('$W_1$ [\textit{J}]', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('PDF', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
%                 
%                 
% figure();
% [fm,xm] = ksdensity(Work(:,2));
% [~, index] = min(abs(xm-W_Sorted_median(1,2)));
% plot(xm,fm,'.-r')
% hold on
% line([W_Sorted_median(1,2) W_Sorted_median(1,2)], [min(fm) max(fm)*(1+0.1) ]);
% % text(W_Sorted_median(1,2)+W_Sorted_median(1,2)*.00003,max(fm),'\rightarrow a_{crit}','FontSize',16);
% hold on
% p = xm(index:end);
% q = fm(index:end);
% H = area(p,q);
% set(H(1),'FaceColor',[1 0.5 0])
% xlabel('$W_2$ [\textit{J}]', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('PDF', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% 
%                 
%                 
% figure();
% [fm,xm] = ksdensity(Work(:,3));
% [~, index] = min(abs(xm-W_Sorted_median(1,3)));
% plot(xm,fm,'.-r')
% hold on
% line([W_Sorted_median(1,3) W_Sorted_median(1,3)], [min(fm) max(fm)*(1+0.1) ]);
% % text(W_Sorted_median(1,3)+W_Sorted_median(1,3)*.00003,max(fm),'\rightarrow _{crit}','FontSize',16);
% hold on
% p = xm(index:end);
% q = fm(index:end);
% H = area(p,q);
% set(H(1),'FaceColor',[1 0.5 0])
% xlabel('$W_3$ [\textit{J}]', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% ylabel('PDF', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');

%% Function: Forman_K_gp
% function [a_all,n,Work,W_Sorted_mean,W_Sorted_median] = ...
%          Forman_K_gp(T,Fmax, Fmin, a0, n0, nf, m,...
%          C_mu, C_cov, N_Samples)

function [a_all,n,W_tot] = ...
         Forman_K_gp(T,Fmax, Fmin, a0, n0, nf, m,...
         C_Samples, N_Samples)
     
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
a_all = zeros(N_Samples, length(n));

% delta_K = zeros(N_Samples, length(n)-1);
% da_dn = zeros(N_Samples, length(n)-1);
a_all(:, 1) = a0(1,1) * ones(N_Samples, 1);

% Work done
W1_k = zeros(N_Samples, length(n)-1,length(Fmax));
W2_k = zeros(N_Samples, length(n)-1,length(Fmax));
W_all = zeros(N_Samples, length(n)-1,length(Fmax));

% % 3 blocks for each mission
% W_tot_1 = zeros(N_Samples, missions);
% W_tot_2 = zeros(N_Samples, missions);
% W_tot_3 = zeros(N_Samples, missions);

nf_new = [0 nf];
nf_cum = cumsum(nf_new);
%     ind = 1;

Fmax = repmat(Fmax,N_Samples,1);
Fmin = repmat(Fmin,N_Samples,1);

C = C_Samples(:);
mu_e = log(1 / sqrt(1 + 0.1^2));
sigma_e = log(1 + 0.1^2);

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
        
        if a_dummy>0.04
            dadn=0;
        else
            % gp model input requirement: force in lbs, crack length in mm
            [Kmax_mean] = predict(GPR_K1, [Fmax(:,k) 1000*a_dummy]);
            [Kmin_mean] = predict(GPR_K1, [Fmin(:,k) 1000*a_dummy]);
            
            deltaK = Kmax_mean - Kmin_mean;
            % deltaK from gp model use MPa*sqrt(m) in the Forman Law
            
            deltaK = deltaK/(1e6);
            
            % determine da/dn value - a here is HALF crack length
            %                 dadn = (C * (deltaK)^m/ ((1 - R(k)) * 67 - deltaK) );% + normrnd(6.505792002037899e-08,1.986698609510128e-08); %+ eps_Forman);% + abs(eps_obs);
            dadn = (C .* (deltaK).^m./ ((1 - R(k)) * 67 - deltaK) )*lognrnd(mu_e,sigma_e);
        end
        
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
    
end

W_tot_1 = sum(W_all(:,:,1),2);
W_tot_2 = sum(W_all(:,:,2),2);
W_tot_3 = sum(W_all(:,:,3),2);

% pass mean of work dones of all samples
W_tot = [mean(W_tot_1); mean(W_tot_2); mean(W_tot_3)];

% % Total amount of work done for each blocks and for different samples
% % for each mission
% kk = 0; jj=0;
% for ii=1:missions
%     jj = jj+1;
%     W_tot_1(:,jj) = sum(W_all(:,:,1+kk),2);
%     W_tot_2(:,jj) = sum(W_all(:,:,2+kk),2);
%     W_tot_3(:,jj) = sum(W_all(:,:,3+kk),2);
%     kk = kk + missions;
% end
% Work = [W_tot_1 W_tot_2 W_tot_3];
% % pass mean of work dones of all samples
% W_tot_mean = mean(Work);
% 
% W_tot_median = median(Work);
% 
% % W_tott=Work
% % W_tot = [W_tott(1,1), W_tott(1,5), W_tott(1,9),W_tott(1,2), W_tott(1,6),...
% %     W_tott(1,10), W_tott(1,3), W_tott(1,7), W_tott(1,11),...
% %     W_tott(1,4), W_tott(1,8), W_tott(1,12)];
% 
% % W_Sorted = [W_tott(1,1), W_tott(1,4), W_tott(1,7),W_tott(1,2), W_tott(1,5),...
% %     W_tott(1,8), W_tott(1,3), W_tott(1,6), W_tott(1,9)];
% W_Sorted_mean = [W_tot_mean(1,1), W_tot_mean(1,2), W_tot_mean(1,3)];
% W_Sorted_median = [W_tot_median(1,1), W_tot_median(1,2), W_tot_median(1,3)];

end