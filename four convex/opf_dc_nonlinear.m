function [v, p, obj, info] = opf_dc_nonlinear(G, pmin, pmax, slack, deltaV, varargin)
% OPF_DC_NONLINEAR — 论文 Model 1（非凸）DC-OPF
% minimize   v'*G*v
% s.t.       p_k / v_k = sum_m G(k,m)*v_m
%            |v-1| <= deltaV,  v(slack)=1
%            pmin <= p <= pmax   (非 slack)
%
% 可选参数：
%   'init', struct('v',v0,'p',p0)   % 热启动
%   'solver','interior-point'|'sqp' % fmincon 算法
%   'display','off'|'iter'

n = size(G,1);
G = (G+G.')/2;               % 数值对称

% 解析可选项
pars.init_v = []; pars.init_p = [];
pars.solver  = 'interior-point';
pars.Display = 'off';
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'init'
            s = varargin{i+1};
            if isfield(s,'v'), pars.init_v = s.v; end
            if isfield(s,'p'), pars.init_p = s.p; end
        case 'solver'
            pars.solver = varargin{i+1};
        case 'display'
            pars.Display = varargin{i+1};
    end
end

% 变量 x = [v; p]，用上下界表达 |v-1|<=deltaV 与 v(slack)=1、p 的箱约束
v0   = 1.0;
lb_v = (1 - deltaV)*ones(n,1);  ub_v = (1 + deltaV)*ones(n,1);
lb_v(slack) = v0;               ub_v(slack) = v0;

lb_p = pmin(:);  ub_p = pmax(:);
lb_p(slack) = -Inf;  ub_p(slack) = Inf;   % slack 功率自由以平衡

% 初值（建议传入 LIN/SOC/SDP/SCA 的结果）
v_init = ones(n,1); v_init(slack) = v0;
p_init = zeros(n,1);
if ~isempty(pars.init_v), v_init = pars.init_v(:); end
if ~isempty(pars.init_p), p_init = pars.init_p(:); end
v_init = min(max(v_init,lb_v),ub_v);
p_init = min(max(p_init,lb_p),ub_p);
x0 = [v_init; p_init];

% 目标与非线性等式约束
objfun  = @(x) x(1:n).'*G*x(1:n);
nonlcon = @(x) deal([], x(n+1:end)./x(1:n) - G*x(1:n));  % c=[], ceq= p./v - G*v

% fmincon 设置
opts = optimoptions('fmincon', ...
    'Algorithm', pars.solver, ...
    'Display',   pars.Display, ...
    'SpecifyObjectiveGradient',    false, ...
    'SpecifyConstraintGradient',   false, ...
    'OptimalityTolerance', 1e-9, ...
    'ConstraintTolerance', 1e-9, ...
    'MaxIterations', 1000, ...
    'StepTolerance', 1e-12);

% 线性/边界约束
A=[]; b=[]; Aeq=[]; beq=[];
lb = [lb_v; lb_p];  ub = [ub_v; ub_p];

[x_opt,fval,exitflag,output] = fmincon(objfun, x0, A,b,Aeq,beq, lb,ub, nonlcon, opts);

v = x_opt(1:n);
p = x_opt(n+1:end);
obj = v.'*G*v;

% 诊断
res = p./v - G*v;
info.exitflag = exitflag;
info.output   = output;
info.res_inf  = norm(res, inf);
info.res_2    = norm(res, 2);
end
