function [v,p,W,obj] = opf_dc_sdp(G,pmin,pmax,slack,deltaV)
n = size(G,1);
cvx_begin sdp quiet
    % cvx_solver sdpt3
    variable W(n,n) symmetric
    variable p(n,1)
    minimize( trace(G*W) )
    subject to
        W == semidefinite(n);          % W ⪰ 0
        W(slack,slack) == 1;           % 恒压
        for k = 1:n
            (1-deltaV)^2 <= W(k,k) <= (1+deltaV)^2;
        end
        diag(G*W) == p;                % p = diag(GW)
        p(2:end) >= pmin(2:end);
        p(2:end) <= pmax(2:end);
cvx_end
% 由最大特征值恢复 v 
[Phi,D] = eig(W); [lambda,ix] = max(diag(D));
v = Phi(:,ix) * sqrt(max(lambda,0));
if v(slack) < 0, v = -v; end
obj = v.'*G*v;
p   = diag(G*W);
end
