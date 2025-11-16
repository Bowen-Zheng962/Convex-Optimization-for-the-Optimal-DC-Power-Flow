function [v,p,W,obj] = opf_dc_soc(G,pmin,pmax,slack,deltaV)
n = size(G,1);
cvx_begin quiet
    % cvx_solver sdpt3  % 或者 Mosek（若已安装）
    variable W(n,n) symmetric
    variable p(n,1)
    minimize( trace(G*W) )
    subject to
        % 电压约束：w_kk = v_k^2, 以及 slack
        W(slack,slack) == 1;
        for k = 1:n
            (1-deltaV)^2 <= W(k,k) <= (1+deltaV)^2;
        end
        % 功率等式：p = diag(G*W)
        diag(G*W) == p;
        % 双曲锥松弛：||[2*w_km ; w_kk - w_mm]||_2 <= w_kk + w_mm
        for k = 1:n
            for m = k:n
                norm([2*W(k,m); W(k,k)-W(m,m)]) <= W(k,k)+W(m,m);
            end
        end
        % 节点功率箱约束（除去并网节点）
        p(2:end) >= pmin(2:end);
        p(2:end) <= pmax(2:end);
cvx_end
% 从对角恢复电压幅值（SOC 用 sqrt(diag(W)) 最稳）
v = sqrt(max(0,diag(W)));
if v(slack) < 0, v = -v; end
obj = v.'*G*v;         % 网损
p   = diag(G*W);       % 一次性返回一致的 p
end
