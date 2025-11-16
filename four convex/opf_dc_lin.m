function [v,p,obj] = opf_dc_lin(G,pmin,pmax,slack,deltaV)
n = size(G,1);
cvx_begin quiet
    variable v(n,1)
    variable p(n,1)
    minimize( quad_form(v,G) )   % v' G v
    subject to
        v(slack) == 1;
        abs(v - 1) <= deltaV;    % |v-1| ≤ Δ
        % 线性化潮流：p_k = sum_m g_km (v_k + v_m - 1)
        for k = 1:n
            p(k) == G(k,:)*( v + v(k)*ones(n,1) - ones(n,1) );
        end
        p(2:end) >= pmin(2:end);
        p(2:end) <= pmax(2:end);
cvx_end
obj = v.'*G*v;
end
