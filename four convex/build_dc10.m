
function [G,pmin,pmax,slack,deltaV,nodenames] = build_dc10()

slack = 1;                 % 节点1为恒压
deltaV = 0.0025;            % |v-1| <= 0.003  (≈0.3%)
n = 10;
nodenames = (1:n)';

% 线路 (i,j) 及 r_ij [pu] —— 按表格顺序
E = [ 1 2;
      2 3;
      2 4;
      4 5;
      2 6;
      6 7;
      7 8;
      7 9;
      3 10 ];
r = [0.0050;
     0.0015;
     0.0020;
     0.0018;
     0.0023;
     0.0017;
     0.0021;
     0.0013;
     0.0015];

% 构造电导矩阵 G
G = zeros(n);
for e = 1:size(E,1)
    i = E(e,1); j = E(e,2); g = 1/r(e);
    G(i,i) = G(i,i) + g;  G(j,j) = G(j,j) + g;
    G(i,j) = G(i,j) - g;  G(j,i) = G(i,j);
end

% 节点功率边界 (pu) —— 表1的 Pmin/Pmax 对应“第二列节点”：
pmin = -inf(n,1);  pmax = +inf(n,1);
pmin(2)  = 0.00;  pmax(2)  = 0.30;
pmin(3)  = -1.80; pmax(3)  = -1.00;
pmin(4)  = -2.30; pmax(4)  = -2.00;
pmin(5)  = 0.00;  pmax(5)  = 3.00;
pmin(6)  = 0.00;  pmax(6)  = 2.00;
pmin(7)  = 0.00;  pmax(7)  = 0.50;
pmin(8)  = 0.00;  pmax(8)  = 1.30;
pmin(9)  = -3.00; pmax(9)  = -1.00;
pmin(10) = 0.00;  pmax(10) = 2.25;
% 节点1(并网)不设上下界
end
