function [Fsol,trials] = BlockLoad(a0,m,n_samples)
%==========================================================================    
%     
%   PARAMETERS:
%    
%   a0: initial crack size
%   a_crit: critical crack size
%   N0: initial # of cyles
%   B: constant value in Paris' law
%   m: material property for Paris law
%   n_samples: # of samples for MCS
%
%   UniAxial Case
%   -------------
%   3 design variables (dv) for each load block (BL)
%   BL: # of block loads
%   # of dv = 3*BL
%
%   Each block has it's own amplitude (A), frequency (w) and time (t)
%==========================================================================    

    format long
    % Load Gaussian Process model for K1 and Work done
    GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
    GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');
    GP = {GPStruc_K1 GPStruc_W};

    %% =================   surrogateopt Inputs  ===========================
    a_crit = 19.337218673249; % mm
    
    % randomly choose Minimum Work done obtained from global opt.
    W_Global = 4.232302423175552e+04;
    t_global = 4219; % [s]
    
    % mission 1 work done assuming 4 missions, dividing into 4 intervals
%     W_total_1 = W_Global*0.30; % 1st mission 30% of total work done
%     t_total_1 = t_global*0.30; % 1st mission 30% of total # of cycles
%     W_total_1 = W_Global*0.20; % 2nd mission 20% of total work done
%     t_total_1 = t_global*0.20; % 2nd mission 20% of total # of cycles
%     W_total_1 = W_Global*0.40; % 3rd mission 40% of total work done
%     t_total_1 = t_global*0.40; % 3rd mission 40% of total # of cycles
    W_total_1 = W_Global*0.10; % 3rd mission 40% of total work done
    t_total_1 = t_global*0.10; % 3rd mission 40% of total # of cycles
    
    % lower and upper bound constraints
    F1_lb=3000; F1_ub=4000;
    F2_lb=5000; F2_ub=5400;
	F3_lb=4000; F3_ub=5000;

%     lb = [F1_lb F2_lb F3_lb 295 395 295];
%     ub = [F1_ub F2_ub F3_ub 305 405 305];
    lb = [F1_lb F2_lb F3_lb 0.2*t_total_1*0.8 0.6*t_total_1*0.8 0.2*t_total_1*0.8];
    ub = [F1_ub F2_ub F3_ub 0.2*t_total_1*1.2 0.6*t_total_1*1.2 0.2*t_total_1*1.2];

    % 1st mission appx 5000 cycles
    
    % fixed frequencies of each block loading
    w1 = 5; w2 = 5; w3 = 5;
    w = [w1 w2 w3];
	
    % Initial value for design variables; amplitude, time
    % [A1 A2 A3 t1 t2 t3]
    x_initial = [3500 5200 4500 0.2*t_total_1 0.6*t_total_1 0.2*t_total_1];

    %% MC Samples
    % % Mean of C, material parameter
    mu_C = 4.5628e-9; cov_C = 0.10;
    
    % icdf of the parameters
    U = lhsdesign(n_samples, 1);
    C_Samples = icdf('normal', U(:, 1), mu_C, mu_C*cov_C)


%     W_total_2 = W_Global*0.30; % 1st mission 30% of total work done
%     W_total_3 = W_Global*0.40; % 1st mission 40% of total work done
%     W_total_4 = W_Global*0.10; % 1st mission 10% of total work done
    %     W_total_1 = rand*(W_Global/4 - 0) + 0
%     W_total_2 = rand*(W_Global/2 - W_Global/4) + W_Global/4;
    
    objectivefun = @(x)objfun(x,w,a0,m,n_samples,C_Samples,GP,W_total_1,t_total_1,a_crit);
                          
     options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
        'InitialPoints',x_initial,'UseParallel',true);%,'MaxFunctionEvaluations',30);
    %,'ObjectiveLimit',1/n_samples-1e-13);

    % Call Optimization
    [Fsol,Fval,ex,out,trials] = surrogateopt(objectivefun,lb,ub,options);
    
    %% %%%%%%%%%%%%%%%%%%     Objective function    %%%%%%%%%%%%%%%%%%%%%%%
    function [Pf] = objfun(x,w,a0,m,n_samples,C_Samples,GP,W_total_1,t_total_1,a_crit) 
        
        % Objective function that needs to be minimized: P[a(Nf) > a_crit]
        % maximum loadings
        Fmax_1=x(1); Fmax_2=x(2); Fmax_3=x(3);
        % min. loadings => R = 0.5
        Fmin_1=0.5*Fmax_1; Fmin_2=0.5*Fmax_2; Fmin_3=0.5*Fmax_3;
        
        % final number of cycles for individual blocks: (freq x time)
        nf_1 = w(1)*x(4); nf_2 = w(2)*x(5); nf_3 = w(3)*x(6);
        
		% collect into arrays
		Fmax = [Fmax_1 Fmax_2 Fmax_3]; Fmin = [Fmin_1 Fmin_2 Fmin_3];
        N_f = [nf_1 nf_2 nf_3];
		
        % Crack propagation
        % NumTimes: # of times the predicted final crack size is greater
        % than the critical crack size
        [W, a_f, NumTimes] = Forman_K_W_gp(GP,Fmax,Fmin,a0,...
                                        m,n_samples,N_f,C_Samples,a_crit);
         
        % Prob. of failure
        Pf = (NumTimes)/n_samples;
        
        % total duration of 3 blocks
        t_total = x(4) + x(5) + x(6);

%         Pf = mean(a_f);

        % penalty function: total work done in that mission is less than
        % the minimum required work done that is randomly chosen from the
        % total work done during the service life, which is calculated
        % using the global optimization
        if ( sum(W) <  W_total_1)
            Pf = 2*Pf;
        end
        % penalty function
        if (t_total < t_total_1) % 20% of the total time 1st mission
           Pf = 2*Pf;
        end

    end


%% %%%%%%%%%%%%%%%%%%%%        Display & Print     %%%%%%%%%%%%%%%%%%%%%%%%
% tmp = load('Work.mat','a_final');
% aF = tmp.('a_final');

% tmp = load('Work.mat','af');
% aF = tmp.('af');

% assignin('base','af',af);
% assignin('base','a_f',a_f)

% % Display and Plotting
fprintf('Amplitude of 1st block: %1.1f [lbs]\n',Fsol(1));
fprintf('Amplitude of 2nd block: %1.1f [lbs]\n',Fsol(2));
fprintf('Amplitude of 3rd block: %1.1f [lbs]\n',Fsol(3));
% fprintf('Frequency of first block: %1.3f [s^-1]\n',Fsol(3));
% fprintf('Frequency of second block: %1.3f [s^-1]\n',Fsol(4));
fprintf('Duration of 1st block loading: %.1f [s]\n',Fsol(4));
fprintf('Duration of 2nd block loading: %.1f [s]\n',Fsol(5));
fprintf('Duration of 3rd block loading: %.1f [s]\n',Fsol(6));
fprintf('# of cycles of the 1st block loading: %.1f\n',round(Fsol(4)*w(1)));
fprintf('# of cycles of the 2nd block loading: %.1f\n',round(Fsol(5)*w(2)));
fprintf('# of cycles of the 3rd block loading: %.1f\n',round(Fsol(6)*w(3)));
fprintf('Total # of cycles: %.1f [s]\n',round(w(1)*Fsol(4)+w(2)*Fsol(5)+w(3)*Fsol(6)));
% fprintf('Probability of failure: %f\n',floor(Fval/a_crit)/n_samples);
fprintf('Objectove function: %f\n',Fval);
% if F>=1
%     fprintf('%d constraint(s) is/(are) not satisfied \n',floor(Fval));
% else
%     fprintf('Probability of failure: %f\n',(Fval-floor(Fval)));
% end


% Plot block loads    
% Define some parameters that define the triangle wave.
elementsPerHalfPeriod1 = Fsol(4)/2; % Number of elements in each rising or falling section.
elementsPerHalfPeriod2 = Fsol(5)/2; % Number of elements in each rising or falling section.
elementsPerHalfPeriod3 = Fsol(6)/2; % Number of elements in each rising or falling section.
% Peak-to-peak amplitude.
fmax1 = Fsol(1); fmax2 = Fsol(2); fmax3 = Fsol(3); % Peak-to-peak amplitude
verticalOffset = -2; % Also acts as a phase shift.
numberOfPeriods1 = w(1); % How many replicates of the triangle you want.
numberOfPeriods2 = w(2); % Basically the frequencies
numberOfPeriods3 = w(3);
% Construct one cycle, up and down.
risingSignal1 = linspace(0.5*fmax1, fmax1, elementsPerHalfPeriod1);
fallingSignal1 = linspace(fmax1, 0.5*fmax1, elementsPerHalfPeriod1);
risingSignal2 = linspace(0.5*fmax2, fmax2, elementsPerHalfPeriod2);
fallingSignal2 = linspace(fmax2, 0.5*fmax2, elementsPerHalfPeriod2);
risingSignal3 = linspace(0.5*fmax3, fmax3, elementsPerHalfPeriod3);
fallingSignal3 = linspace(fmax3, 0.5*fmax3, elementsPerHalfPeriod3);
% Combine rising and falling sections into one single triangle.
firstCycle = [risingSignal1, fallingSignal1(2:end)] + verticalOffset; 
secondCycle = [risingSignal2, fallingSignal2(2:end)] + verticalOffset;
thirdCycle = [risingSignal3, fallingSignal3(2:end)] + verticalOffset; 
% Now replicate this cycle several (numberOfPeriods) times.
waveform1 = repmat(firstCycle, [1 round(numberOfPeriods1)]);
x1 = 0 : length(waveform1)-1;
waveform2 = repmat(secondCycle, [1 round(numberOfPeriods2)]);
x2 = x1(end) : x1(end)+length(waveform2)-1;
waveform3 = repmat(thirdCycle, [1 round(numberOfPeriods3)]);
x3 = x2(end) : x2(end)+length(waveform3)-1;
% Now plot the triangle wave.
figure();
n0=27500; % initial # of cycles
plot(n0+x1,waveform1,'b-',n0+x2,waveform2,'r-',n0+x3,waveform3,'g-','LineWidth',2);
grid on;
title('Load Spectrum', 'FontName', 'Times New Roman', ...
                    'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')
xlabel('Cycles', 'FontName', 'Times New Roman', ...
                    'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')
ylabel('Load [\textit{lbf}]', 'FontName', 'Times New Roman', ...
                    'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')


    
    
    
% %% Plot PDF of crack size
% figure();
% a_F = a_F * 1e3; % in [mm]
% a_crit = a_crit*1e3; % in [mm]
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
% xlim([a0*1e3 a_crit+2e-1])

    

end