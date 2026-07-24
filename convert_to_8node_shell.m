function [nodes_ext, quads8] = convert_to_8node_shell(nodes, quads)
% 将四节点单元拓展为八节点壳单元（共享边中点只生成一次）
% 输入：
%   nodes: 原始坐标 N x 3
%   quads: 原始四边形单元 M x 4 [n1 n2 n3 n4]
% 输出：
%   nodes_ext: 扩展节点坐标（含中点），(N + ?) x 3
%   quads8: M x 8, 每行为 [n1, n2, n3, n4, m12, m23, m34, m41]

nodes_ext = nodes;
node_count = size(nodes, 1);
n_quads = size(quads, 1);
quads8 = zeros(n_quads, 8);

% 用来记录已有的边中点节点：key = 'minID_maxID'，value = 新节点编号
edge_mid_map = containers.Map();

for i = 1:n_quads
    n = quads(i, 1:4); % 节点顺序 n1, n2, n3, n4
    mids = zeros(1, 4); % 保存四条边的中点节点编号
    
    % 定义边对：1-2, 2-3, 3-4, 4-1
    edge_pairs = [1 2; 2 3; 3 4; 4 1];
    
    for k = 1:4
        a = n(edge_pairs(k,1));
        b = n(edge_pairs(k,2));
        key = sprintf('%d_%d', min(a,b), max(a,b)); % 唯一边标识
        
        if isKey(edge_mid_map, key)
            mids(k) = edge_mid_map(key); % 已存在，直接复用
        else
            % 计算中点，添加新节点
            p_mid = (nodes(a,:) + nodes(b,:)) / 2;
            nodes_ext = [nodes_ext; p_mid];
            node_count = node_count + 1;
            mids(k) = node_count;
            edge_mid_map(key) = node_count;
        end
    end
    
    % 组装八节点单元行：[n1, n2, n3, n4, m12, m23, m34, m41]
    quads8(i, :) = [n(1), n(2), n(3), n(4), mids(1), mids(2), mids(3), mids(4)];
end

end