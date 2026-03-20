%% 100-LINE MATLAB SOURCE CODES of BALANCING COMPOSITE MOTION OPTIMIZATION (BCMO)
%% AUTHORS: Thang Le-Duc, Quoc-Hung Nguyen, H. Nguyen-Xuan
%% Start main program
clear all, close all, clc
%% Define input parameters
total_time = 10;    % Total times for solving problem
Function_name='F1'; % Name of the test function 
% Call benchmark functions collected from papers "Evolutionary Programming Made
% Faster" and "Grey Wold Optimizer"

[lb,ub,dim,Objf]=Get_Functions_details(Function_name); 

LB=repmat(lb,1,dim); UB=repmat(ub,1,dim);

NP = 50;                   % Population size
MaxGen = 1000;             % Maximum number of iterations

%% BCMO algorithm for global optimization problem
tpad = cputime;               % Start to calculate the elapsed time
for time=1:total_time
    [optimalSol_BCMO,BCMO_Converg{time}] = BCMO(Objf,NP,MaxGen,dim,LB,UB);%Call BCMO
    resultBCMO(time,:) = optimalSol_BCMO;  
end

%% Save the optimal result
totaltime = (cputime-tpad)      % Print the elapsed time
[bestValueBCMO,id] = min(resultBCMO(:,dim+1)); % Print the best objective function value

worst = max(resultBCMO(:,dim+1))      % Print the worst objective function value
Ave_ = mean(resultBCMO(:,dim+1))      % Print the average objective function value
Std_ = std(resultBCMO(:,dim+1))       % Print the standard deviation of objective function value

%save resultBCMO.mat resultBCMO Std_;  % Save all information about optimal result
%save time.mat totaltime;            % Save the elapsed time

display(['Best solution: ', num2str(resultBCMO(id,1:dim))]);
display(['Best optimal value: ', num2str(bestValueBCMO)]);

figure('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1,x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(BCMO_Converg{id},'Color','b')
title('Objective space')
xlabel('Iteration');
ylabel('Fitness function');
axis tight
grid on
box on
legend('Solution covergence')