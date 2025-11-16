clear; clc;
[G,pmin,pmax,slack,deltaV,nodenames] = build_dc10();

% —— SOC
[vsoc,psoc,~,obj_soc] = opf_dc_soc(G,pmin,pmax,slack,deltaV);

% —— SDP
[vsdp,psdp,~,obj_sdp] = opf_dc_sdp(G,pmin,pmax,slack,deltaV);

% —— 线性化
[vlin,plin,obj_lin]   = opf_dc_lin(G,pmin,pmax,slack,deltaV);

% —— SCA (新增)
[vsc, psc, obj_sca, iter_sca] = opf_dc_sca(G, pmin, pmax, slack, deltaV, 'verbose', true);

% === Nonlinear (论文 Model 1 非凸原始OPF) — 用线性化做热启动 ===
init_nl.v = vlin;
init_nl.p = plin;
[vnl, pnl, obj_nl, infonl] = opf_dc_nonlinear(G,pmin,pmax,slack,deltaV, ...
    'init', init_nl, 'solver','interior-point', 'display','off');
fprintf('[NONLIN] obj=%.6f, ||res||_inf=%.2e, exitflag=%d\n', ...
    obj_nl, infonl.res_inf, infonl.exitflag);

% 结果汇总
T = table(nodenames, vsoc, vsdp, vlin, vsc, vnl, ...
                     psoc, psdp, plin, psc, pnl, ...
    'VariableNames', {'Node','V_SOC','V_SDP','V_LIN','V_SCA','V_NL', ...
                             'P_SOC','P_SDP','P_LIN','P_SCA','P_NL'});
disp('=== Results (Voltages & Powers) ===');
disp(T);

fprintf('\nLoss (objective) — SOC: %.6f, SDP: %.6f, LIN: %.6f, SCA: %.6f, NONLIN: %.6f\n', ...
    obj_soc, obj_sdp, obj_lin, obj_sca, obj_nl);
fprintf('SCA 迭代次数: %d\n', iter_sca);
