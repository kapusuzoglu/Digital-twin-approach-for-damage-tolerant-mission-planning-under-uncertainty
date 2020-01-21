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
%             Work
% -----------------------------
% arrange all the data
X_W = zeros(nl*nc, 2); 
% 1st col: load in lbs, 2nd col: crack length in mm
Y_W = zeros(nl*nc, 1);
% col: Work done in J

% plate geometry
lx = 0.0254 * 16; % m, 1 in = 0.0254 m
ly = 0.0254 * 6; % m
t = 0.0254 * 0.032; % m

% structured mesh on the lefr edge, QUADRATIC elements (3 nodes on edge)...
% element size is ly/10000
% size_elem = ly/10000;
size_elem = 0.0005;

% get the list of coordinates of node on the left edge... 
% where the traction is applied, sort and get indecies
list_node = dlmread('outNodes.txt');
y_coord = list_node(:, 3);
[y_coord_sorted, idx] = sort(y_coord);
num_nodes = length(y_coord);

% a list of areas of the traction applied on at each node
Area_traction = 0.5 * t*size_elem * [0.5; ones(num_nodes - 2, 1);0.5];

k = 1;
for i = 1:nl % load = loadarray(i) in lbs
    Force_total = 4.4482 * loadarray(i); % in Newtons
    Traction = Force_total/(3 * 0.0254 * t); % in Pa
    for j = 1:nc % crack length = clarray(j) in mm
        filename = ['F_' num2str(loadarray(i)) '_CL_' ...
                    num2str(clarray(j)) '_U1HistoryOutput.txt'];
        disp = dlmread(filename); % displacement of each node
        disp_sorted = - disp(idx); % sort according to sorted node list
        % orginal data: negative value, meaning u1 to the left, in meters
        Y_W(k) = sum(Traction * Area_traction .* disp_sorted);
        X_W(k, 1) = loadarray(i);
        X_W(k, 2) = clarray(j);
        k = k + 1;
    end
end


X_W_train = X_W;
Y_W_train = Y_W;

save('TrainData_W.mat',...
     'X_W_train', 'Y_W_train')

% initialize the kernel parameters
sigma0_W = std(Y_W_train);
sigmaF0_W = sigma0_W;
d = size(X_W_train,2);
sigmaM0_W = 10*ones(d,1);

gprMdl_W = fitrgp(X_W_train,Y_W_train,'Basis','constant',...
                   'FitMethod','exact','PredictMethod','exact',...
                   'KernelFunction','ardsquaredexponential',...
                   'KernelParameters',[sigmaM0_W;sigmaF0_W],...
                   'Sigma',sigma0_W,'Standardize',1);

% save the GP model
save('gprMdl_W.mat', 'gprMdl_W')

% Model parameters
% Trend Function
Beta_W= gprMdl_W.Beta;
% The kernel parameters. 
sigmaM_W = gprMdl_W.KernelInformation.KernelParameters(1:end-1,1);
sigmaF_W = gprMdl_W.KernelInformation.KernelParameters(end);
% Estimated noise standard deviation (sigma_n)
sigma_W  = gprMdl_W.Sigma; 
% save the model parameters
save('Param_gprMdl_W.mat',...
     'Beta_W', 'sigmaM_W', 'sigmaF_W', 'sigma_W')

% test the GP model
% W_pred = predict(gprMdl_W,X_W_test);

%%
% -----------------------------
%   Plot GP model test result
% -----------------------------
% figure()
% 
% plot(1:length(W_pred), Y_W_test, 'o', 1:length(W_pred), W_pred, 'x')
% legend('True Value', 'GP Prediction')
% xlabel('Test Point')
% ylabel('$Work\ (J)\ Done$',...
%        'interpreter','latex','fontsize', 14)
% title('Test Result - GP Model for Work Done','fontsize', 16)