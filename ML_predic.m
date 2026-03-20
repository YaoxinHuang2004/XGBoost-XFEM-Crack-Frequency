%% ========================================================================
%  BCMO-XGBoost Batch Prediction
% =========================================================================
clc; clear; close all;

%% ==================== 0. Retrieve and print hardware (CPU) information ====================
try
    if ispc
        [~, cpu_info] = system('wmic cpu get name');
        cpu_info = strtrim(regexprep(strrep(cpu_info, 'Name', ''), '\s+', ' '));
    elseif ismac
        [~, cpu_info] = system('sysctl -n machdep.cpu.brand_string');
        cpu_info = strtrim(cpu_info);
    else
        cpu_info = 'Unknown CPU';
    end
catch
    cpu_info = 'Unable to retrieve CPU info';
end
fprintf('Test platform hardware configuration (CPU): %s\n', cpu_info);
fprintf('--------------------------------------------------\n');

%% ==================== 1. Reading datasets and loading models ====================
filename = 'New_3000_LHS_Dataset.csv';
fprintf('Reading file: %s ...\n', filename);
data = readtable(filename);
num_samples = height(data);

fprintf('Loading BCMO-XGBoost Model...\n');
load('Best_BCMO_XGBoost_Model.mat', 'mdl');

%% ==================== 2. Construct the predictive feature matrix ====================
X_batch = [data.Type_of_crack, ...
           data.Crack_length, ...
           data.a, ...
           data.b, ...
           data.t_b, ...
           data.Type_of_BCs, ...
           data.kx, ...
           data.ky];

fprintf('\n=== ML Batch processing of the model has begun. There are %d data sets in total. ===\n', num_samples);

%% ==================== 3. Batch forecasting ====================
tic; 
% Predict 3,000 results at the touch of a button
Omega_bar_pred = predict(mdl, X_batch);

batch_time = toc; 
fprintf('--------------------------------------------------\n');
fprintf('Batch prediction complete %d natural frequencies have been successfully calculated!\n', num_samples);
fprintf('Total time taken = %.6f s !!!\n', batch_time);
fprintf('Average duration per session = %.8f s\n', batch_time / num_samples);
fprintf('--------------------------------------------------\n');

%% ==================== 4. Save results ====================
% Add the forecast results to the table as a new column
data.Omega_bar_Pred = Omega_bar_pred;

% Write the results to a new Excel file
output_filename = 'ML_Predicted_3000_LHS_Dataset.xlsx';
writetable(data, output_filename);
fprintf('\n*** All calculations are complete! The forecast results have been saved to the last column of %s ***\n', output_filename);