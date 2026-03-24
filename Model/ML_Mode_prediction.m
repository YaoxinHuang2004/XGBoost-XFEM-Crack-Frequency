clear; clc;
% 1. 直接读取 mat 文件中的所有模型和变量
load('v5_model.mat');

fprintf('\n--- STEP 4: NN Prediction Visualization (Mode %d) ---\n\n', target_mode);
[Xp, Yp_g] = meshgrid(linspace(0,1,ngx), linspace(0,1,ngy));

% 8 个测试用例
test_cases = {
    [0, 0.3, 1, 0.02, 0, 1, 1, 1],  'Edge, SSSS, k=1, lc/a=0.3';
    [0, 0.5, 1, 0.02, 1, 1, 1, 1],  'Edge, CCCC, k=1, lc/a=0.5';
    [1, 0.3, 1, 0.02, 0, 1, 1, 1],  'Center, SSSS, k=1, lc/a=0.3';
    [1, 0.5, 1, 0.02, 1, 1, 1, 1],  'Center, CCCC, k=1, lc/a=0.5';
    [0, 0.5, 1, 0.02, 0, 5, 5, 1],  'Edge, SSSS, k=5, lc/a=0.5';
    [0, 0.7, 1, 0.02, 1, 10,10,1],  'Edge, CCCC, k=10, lc/a=0.7';
    [1, 0.5, 2, 0.05, 0, 5, 5, 1],  'Center, SSSS, k=5, rect';
    [0, 0.5, 2, 0.05, 1, 10,10,1],  'Edge, CCCC, k=10, rect';
};

% ====== MAIN FIGURE: 2x4 布局, 仅显示模型预测 3D 曲面 ======
fig = figure('Position',[10,10,1600,800],'Color','w');
for tc = 1:8
    p8 = test_cases{tc,1};
    lbl = test_cases{tc,2};
    
    % ---------------------------------------------------------
    % 用矩阵运算直接替代 nn_predict_v5 函数：
    % 1. 输入特征归一化
    X_norm = (p8 - Xmin) ./ Xrng;
    
    % 2. 神经网络前向传播 (注意 net 需要列向量输入)
    Ys_pred_norm_col = net(X_norm'); 
    Ys_pred_norm = Ys_pred_norm_col'; % 转回行向量
    
    % 3. 反归一化 PCA 得分 (Z-score 逆变换)
    Ys_pred = Ys_pred_norm .* Ys_std + Ys_mu;
    
    % 4. PCA 逆变换恢复完整节点位移 (Y = Score * V' + Mu)
    Y_pred_1D = Ys_pred * V_pca' + Y_mu;
    
    % 5. 重塑为二维物理网格
    Wnn = reshape(Y_pred_1D, [ngx, ngy]);
    % ---------------------------------------------------------

    % 计算裂纹线坐标用于绘图
    lc_v = p8(2);
    if p8(1) == 0
        cx = linspace(0, lc_v, 50);
    else
        cx = linspace(0.5 - lc_v/2, 0.5 + lc_v/2, 50);
    end
    cy_v = 0.5 * ones(size(cx));
    
    % 绘制 3D 预测结果
    subplot(2, 4, tc);
    surf(Xp, Yp_g, Wnn, 'EdgeColor', 'none', 'FaceAlpha', 0.92); hold on;
    % 将裂纹线投影到曲面上
    czn = interp2(Xp, Yp_g, Wnn, cx, cy_v, 'linear', 0);
    plot3(cx, cy_v, czn, 'r-', 'LineWidth', 2.5); hold off;
    
    colormap(jet); colorbar('FontSize', 5);
    title(sprintf('BCMO-XGBR Prediction\n%s', lbl), 'FontSize', 9);
    view([-30, 30]); axis tight; grid on; set(gca, 'FontSize', 6);
end
sgtitle('Mode 1 Shape Prediction: BCMO-XGBR', 'FontSize', 14);


% ====== CONTOUR COMPARISON: 2x4 布局, 仅显示等高线 ======
fig2 = figure('Position',[10,10,1600,800],'Color','w');
for tc = 1:8
    p8 = test_cases{tc,1}; 
    lbl = test_cases{tc,2};
    
    % --- 同样的预测逻辑再跑一遍画 2D 图 ---
    X_norm = (p8 - Xmin) ./ Xrng;
    Ys_pred_norm = (net(X_norm'))'; 
    Ys_pred = Ys_pred_norm .* Ys_std + Ys_mu;
    Y_pred_1D = Ys_pred * V_pca' + Y_mu;
    Wnn = reshape(Y_pred_1D, [ngx, ngy]);
    % -------------------------------------
    
    lc_v = p8(2);
    if p8(1) == 0
        xs = 0; xe = lc_v;
    else
        xs = 0.5 - lc_v/2; xe = 0.5 + lc_v/2;
    end
    
    subplot(2, 4, tc);
    contourf(Xp, Yp_g, Wnn, 15, 'LineColor', 'none'); hold on;
    plot([xs, xe], [0.5, 0.5], 'r-', 'LineWidth', 2); hold off;
    
    colormap(jet); axis equal tight; set(gca, 'FontSize', 5);
    title(sprintf('BCMO-XGBR\n%s', lbl), 'FontSize', 8);
end
sgtitle('Mode 1 Contour Prediction: BCMO-XGBR', 'FontSize', 14);

fprintf('\n============================================================\n');
fprintf('  Done! All plots generated successfully.\n');
fprintf('============================================================\n');