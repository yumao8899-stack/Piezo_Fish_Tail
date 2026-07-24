function [nodes, elements, index] = read_comsol8(filename)
% READ_BDF_INTERFACE 读取 COMSOL 导出的 BDF 网格文件
% 输入:
%   filename - .bdf 文件路径
% 输出:
%   nodes    - N x 3 矩阵，节点坐标 [x, y, z] (z通常为0)
%   elements - M x 8 矩阵，四边形单元的节点编号 (基于 nodes 的行索引，用于 patch 绘图)
%   index    - K x 1 元胞数组 (Cell Array)，K 为最大域编号 (例如 6)
%              index{i} 包含所有属于域 i 的原始单元序号 (Element ID)

    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    raw_nodes = [];     % [ID, X, Y, Z]
    raw_elems = [];     % [EID, PID, N1...N8]
    
    % 正则表达式：修复粘连负号 (例如 1.234-5.678)
    regex_fix_sign = '(?<=[0-9])-(?=[0-9])'; 

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        
        % 跳过空行或注释
        if isempty(line) || line(1) == '$'
            continue; 
        end
        
        % --- 读取节点 (GRID) ---
        if startsWith(line, 'GRID')
            line = regexprep(line, regex_fix_sign, ' -');
            data = sscanf(line(5:end), '%f');
            if length(data) >= 4
                raw_nodes = [raw_nodes; data(1), data(2:4)'];
            end
            
        % --- 读取单元 (CQUAD8) ---
        elseif startsWith(line, 'CQUAD8')
            % 处理续行
            current_str = line;
            while contains(current_str, '+')
                nextLine = fgetl(fid);
                if ~ischar(nextLine), break; end
                nextLine = strtrim(nextLine);
                current_str = strrep(current_str, '+CONT', '');
                current_str = strrep(current_str, '+', ''); 
                nextLine = strrep(nextLine, '+CONT', '');
                nextLine = strrep(nextLine, '+', '');
                current_str = [current_str ' ' nextLine];
            end
            
            current_str = regexprep(current_str, regex_fix_sign, ' -');
            data = sscanf(current_str(7:end), '%f');
            
            % EID(1), PID(2), Nodes(3:10)
            if length(data) >= 10
                raw_elems = [raw_elems; data(1), data(2), data(3:10)'];
            end
        end
    end
    fclose(fid);

    if isempty(raw_nodes) || isempty(raw_elems)
        error('未读取到有效数据');
    end

    % --- 1. 处理 Nodes ---
    % 按 ID 排序 (可选，但推荐)
    [~, sortIdx] = sort(raw_nodes(:,1));
    sorted_nodes = raw_nodes(sortIdx, :);
    nodes = sorted_nodes(:, 2:4); % 输出坐标
    
    % 建立 ID 映射表: 原始ID -> MATLAB行索引
    max_node_id = max(raw_nodes(:,1));
    id_map = zeros(max_node_id, 1);
    id_map(raw_nodes(:,1)) = 1:size(raw_nodes, 1);
    
    % --- 2. 处理 Elements ---
    % 将单元连接关系中的原始节点ID 替换为 行索引
    elements = id_map(raw_elems(:, 3:10));
    
    % --- 3. 处理 Index (按域分组) ---
    all_eids = raw_elems(:, 1); % 第一列：单元序号
    all_pids = raw_elems(:, 2); % 第二列：域编号 (1~6)
    
    max_pid = max(all_pids);    % 应该是 6
    index = cell(max_pid, 1);   % 创建 6x1 的元胞数组
    
    for i = 1:max_pid
        % 找到所有 PID 等于 i 的行，提取对应的 EID
        current_eids = all_eids(all_pids == i);
        % 如果该域有单元，则存入；没有则为空
        if ~isempty(current_eids)
            index{i} = current_eids;
        end
    end

end