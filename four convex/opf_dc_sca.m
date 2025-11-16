function [v, p, obj, iter] = opf_dc_sca(G, pmin, pmax, slack, deltaV, varargin)
% SCA: Successive Convex Approximation for DC OPF
% 输入：
%   G: 导纳矩阵 (n x n)
%   pmin, pmax: 功率上下限 (n x 1)
%   slack: 平衡节点索引
%   deltaV: 电压偏差上限 |v-1| <= deltaV
%   可选参数：'tol', 'max_iter', 'verbose'
%
% 输出：
%   v: 电压
%   p: 功率
%   obj: 网损 v'*G*v
%   iter: 迭代次数

% 默认参数
tol = 1e-6;
max_iter = 50;
verbose = true;

% 解析可选参数
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'tol'
            tol = varargin{i+1};
        case 'max_iter'
            max_iter = varargin{i+1};
        case 'verbose'
            verbose = varargin{i+1};
    end
end

n = size(G,1);

% 初始化：平启动
v = ones(n,1);  % pu
converged = false;

if verbose
    fprintf('SCA 迭代开始...\n');
    fprintf('%3s %12s %12s\n', 'Iter', 'Objective', '||Δv||');
end

for iter = 1:max_iter
    v_old = v;
    
    % 计算当前电流 I_k = sum_m g_km * v_m
    I = G * v;  % 节点注入电流
    
    % 开始 CVX 优化
    cvx_begin quiet
        variable v_new(n)
        variable p_new(n)
        
        minimize( quad_form(v_new, G) )  % v_new' * G * v_new
        
        subject to
            % Slack 节点电压固定
            v_new(slack) == 1;
            
            % 电压范围约束
            abs(v_new - 1) <= deltaV;
            
            % SCA 线性化潮流约束：p_k ≈ v_k^{old} * I_k + v_k * I_k^{old} - v_k^{old} I_k^{old}
            % 其中 I_k = sum_m g_km * v_m_new
            I_new = G * v_new;  % 新电流
            for k = 1:n
                p_new(k) == v_old(k) * I_new(k) + v_new(k) * I(k) - v_old(k) * I(k);
            end
            
            % 功率上下限（除平衡节点外）
            p_new(slack) == [];  % 平衡节点功率自由
            p_new([1:slack-1, slack+1:n]) >= pmin([1:slack-1, slack+1:n]);
            p_new([1:slack-1, slack+1:n]) <= pmax([1:slack-1, slack+1:n]);
    cvx_end
    
    if ~exist('cvx_status','var') || strcmp(cvx_status,'Infeasible') || strcmp(cvx_status,'Failed')
        error('SCA: 内层凸问题求解失败，状态 = %s', cvx_status);
    end
    
    v = v_new;
    p = p_new;
    
    % 计算变化量
    dv = norm(v - v_old, 'inf');
    obj = v' * G * v;
    
    if verbose
        fprintf('%3d %12.6f %12.2e\n', iter, obj, dv);
    end
    
    if dv < tol
        converged = true;
        break;
    end
end

if ~converged && verbose
    fprintf('SCA 警告：未在 %d 次迭代内收敛。\n', max_iter);
end

end