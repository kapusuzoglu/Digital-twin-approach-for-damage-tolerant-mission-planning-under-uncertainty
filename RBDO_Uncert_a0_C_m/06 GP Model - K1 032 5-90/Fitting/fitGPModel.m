% Get the Current Directory
Directory = pwd;
% change to current directory
cd(Directory);
% load array
loadarray = 1000:1000:8000;
nl = length(loadarray);
% crack length array
clarray = 5:5:90;
nc = length(clarray);

%%
% -----------------------------
%              K1
% -----------------------------
% arrange all the data
X_K1 = zeros(nl*nc, 2); 
% 1st col: load in lbs, 2nd col: crack length in mm
Y_K1 = zeros(nl*nc, 1);
% col: K1 in ( Pa * sqrt(m) )

k = 1;
for i = 1:nl % load = loadarray(i) in lbs
    for j = 1:nc % crack length = clarray(j) in mm
        filename = ['F_' num2str(loadarray(i)) '_CL_' ...
                    num2str(clarray(j)) '_K1.txt'];
        temp = dlmread(filename);
        Y_K1(k) = temp(1, 3);
        X_K1(k, 1) = loadarray(i);
        X_K1(k, 2) = clarray(j);
        k = k + 1;
    end
end

% randomly select 90% to be training set for GP model
num_train = nl*nc;

% Use the first 90% to train the GP model

X_K1_train = X_K1;
Y_K1_train = Y_K1;

% save the training and testing data
save('TrainData_K1.mat', 'X_K1_train', 'Y_K1_train')

% initialize the kernel parameters
sigma0_K1 = std(Y_K1_train);
sigmaF0_K1 = sigma0_K1;
d = size(X_K1_train,2);
sigmaM0_K1 = 10*ones(d,1);

% fitting using fitrgp
gprMdl_K1 = fitrgp(X_K1_train,Y_K1_train,'Basis','constant',...
                   'FitMethod','exact','PredictMethod','exact',...
                   'KernelFunction','ardsquaredexponential',...
                   'KernelParameters',[sigmaM0_K1;sigmaF0_K1],...
                   'Sigma',sigma0_K1,'Standardize',1);
% save the GP model
save('gprMdl_K1.mat', 'gprMdl_K1')

% Model parameters
% Trend Function
Beta_K1= gprMdl_K1.Beta;
% The kernel parameters. 
sigmaM_K1 = gprMdl_K1.KernelInformation.KernelParameters(1:end-1,1);
sigmaF_K1 = gprMdl_K1.KernelInformation.KernelParameters(end);
% Estimated noise standard deviation (sigma_n)
sigma_K1  = gprMdl_K1.Sigma; 
% save the model parameters
save('Param_gprMdl_K1.mat',...
     'Beta_K1', 'sigmaM_K1', 'sigmaF_K1', 'sigma_K1')

% % test the GP model
% K1_pred = predict(gprMdl_K1,X_K1_test);

%%
% % -----------------------------
% %   Plot GP model test result
% % -----------------------------
% figure()
% 
% plot(1:length(K1_pred), Y_K1_test, 'o', 1:length(K1_pred), K1_pred, 'x')
% legend('True Value', 'GP Prediction')
% xlabel('Test Point')
% ylabel('$Stress\ Intensity\ Factor\ K_1\ (MPa \sqrt{m}$)',...
%        'interpreter','latex','fontsize', 14)
% title('Test Result - GP Model for K_1','fontsize', 16)