function [W_tot, a_f, NumTimes] = Forman_K_W_gp(T,Fmax,Fmin,a0,...
                                        m,N_Samples,N_f,C_Samples,a_crit)
% function [NumTimes, W_tot, a_f] = Forman_K_W_gp(T,Fmax,Fmin,a0,a_crit,...
%                                         m,N_Samples,N_f,C_Samples)
     

% ------ fitrgp ------
% T = load('..\GPModel_K1_and_Work','gprMdl_K1');
GPR_K1 = T{1}.('gprMdl_K1');
GPR_W = T{2}.('gprMdl_W');
% ------ fitrgp ------

% ------ Physical Parameters ------
% Fmax_i, Fmin_i, where i refers to the block loading sequence in [lbs]

R_1 = Fmin(1)/Fmax(1); % Strss Ratio
R_2 = Fmin(2)/Fmax(2); % Strss Ratio
R_3 = Fmin(3)/Fmax(3); % Strss Ratio
R = [R_1 R_2 R_3];

% change in # of cyclec in an arrray
N_all = [N_f(1) N_f(2) N_f(3)];

%% Place holders 
% for estimated crack size
% each row is a trace of crack growth corresponding to a set of MC Samples
a_all = zeros(N_Samples,round(max(N_all))+1);
% convert a from meters to micrometers
a_all(:, 1) = a0 * ones(N_Samples, 1);

% Work done
W1_k = zeros(N_Samples,round(max(N_all)),length(Fmax));
W2_k = zeros(N_Samples,round(max(N_all)),length(Fmax));
W_all = zeros(N_Samples,round(max(N_all)),length(Fmax));

% save final crack size for each sample
a_f = zeros(N_Samples,1);

Nall = zeros(1,length(Fmax));
for k = 1:length(Fmax)
    if (round(N_all(k))-N_all(k))>0
        Nall(k) = N_all(k)+0.5+1e-12;
    else
        Nall(k) = N_all(k)-0.5+1e-12;
    end
end

%% =============    Monte Carlo Simulation (MCS)   ================
% Monte Carlo Simulation for uncertain parameters:
% NumTimes = 0;

% MCS
ind = 1;

a_all(:,ind) = a0;

% MCS
C = C_Samples(:);
    
Fmax = repmat(Fmax,N_Samples,1);
Fmin = repmat(Fmin,N_Samples,1);

% Do for each block loading
for k = 1:length(R)
    
    % Calculate crack growth within each blocks
    % pass the final crack size of the current block
    a_all(:,1) = a_all(:,ind);
    
    i = 0;
    while i < Nall(k)
        i=i+1;
        % pass the final crack size of the previous block as
        % initial crack size of the next block
        a_dummy = a_all(:,i);
            
        % gp model input requirement: force in lbs, crack length in mm
        [Kmax] = predict(GPR_K1, [Fmax(:,k) a_dummy]);
        [Kmin] = predict(GPR_K1, [Fmin(:,k) a_dummy]);
        
        % Delta_K
        deltaK = Kmax - Kmin;
        
        % deltaK from gp model use MPa*sqrt(m) in the Forman Law
        deltaK = deltaK/(1e6);
        
        % determine da/dn value - a here is HALF crack length
        % and convert it to mm
        dadn = 1e3*( C .* (deltaK).^m./ ((1 - R(k)) * 67 - deltaK));
        
        if (sum(a_dummy>25))>=1
            idx = a_dummy>25;
            dadn(idx,1)=0;
        end
        
        % Predict Crack size in the next cycle
        a_all(:, i + 1) = a_dummy + 2 * dadn;
        
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
    ind = i+1;
end

% count the times af > a_crit
NumTimes =sum(a_all(:,ind) > a_crit);

a_f(:,1) = a_all(:,ind);

% Total amount of work done for each blocks and for different samples
W_tot_1 = sum(W_all(:,:,1),2);
W_tot_2 = sum(W_all(:,:,2),2);
W_tot_3 = sum(W_all(:,:,3),2);

% pass mean of work dones of all samples
W_tot = [mean(W_tot_1); mean(W_tot_2); mean(W_tot_3)];



end
