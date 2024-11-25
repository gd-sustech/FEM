% 一维弹性杆有限元分析
clc; clear;

%% 参数设置
123
L = 10;          % 杆的总长度 (m)
E = 200e9;       % 弹性模量 (Pa)
A = 0.01;        % 截面积 (m^2)
F = 1000;        % 右端受力 (N)
n_elem = 2;      % 单元数量
n_nodes = n_elem + 1;  % 节点数量
node_coords = linspace(0, L, n_nodes); % 节点坐标

%% 单元划分与刚度矩阵组装
K = zeros(n_nodes, n_nodes); % 初始化全局刚度矩阵

% 每个单元长度
le = L / n_elem;

% 单元刚度矩阵 (1D 梁公式: k = EA/le * [1 -1; -1 1])
k_elem = (E * A / le) * [1 -1; -1 1];

for i = 1:n_elem
    K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + k_elem;
end

%% 边界条件与载荷向量
F_vec = zeros(n_nodes, 1); % 初始化载荷向量
F_vec(end) = F;            % 在右端节点施加力

% 施加位移边界条件 (左端固定: u(1) = 0)
K_reduced = K(2:end, 2:end);
F_reduced = F_vec(2:end);

%% 求解方程组
u_reduced = K_reduced \ F_reduced; % 位移向量
u = [0; u_reduced];                % 添加固定边界位移

%% 输出结果
disp('节点坐标 (m):');
disp(node_coords');
disp('节点位移 (m):');
disp(u);
