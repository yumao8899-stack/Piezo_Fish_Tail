function [quadNodes, quadElements,  index] = read_comsol4(filename)
% 读取COMSOL网格文件，输出四边形单元的节点坐标和编号
% quadNodes: N x 3 节点坐标，第三列全为0
% quadElements: M x 4 四边形单元的节点编号（Matlab从1开始）

fid = fopen(filename, 'r');
if fid < 0
    error('无法打开文件');
end

% 读取节点数
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, '# number of mesh vertices')
        numNodes = sscanf(tline, '%d');
        break;
    end
end

% 跳到节点坐标部分
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, '# Mesh vertex coordinates')
        break;
    end
end

% 读取节点坐标
quadNodes = zeros(numNodes, 3);
for i = 1:numNodes
    tline = fgetl(fid);
    vals = sscanf(tline, '%f');
    quadNodes(i,1:2) = vals(1:2)';
end

% 跳到四边形单元部分
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, 'quad # type name')
        break;
    end
end

% 读取单元节点数
while ~feof(fid)
    tline = fgetl(fid);
    if ~isempty(sscanf(tline, '%d'))
        numNodesPerElem = sscanf(tline, '%d');
        break;
    end
end

% 读取单元数
while ~feof(fid)
    tline = fgetl(fid);
    if ~isempty(sscanf(tline, '%d'))
        numElems = sscanf(tline, '%d');
        break;
    end
end

% 跳到 # Elements
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, '# Elements')
        break;
    end
end

% 读取单元
quadElements = zeros(numElems, numNodesPerElem);
for i = 1:numElems
    tline = fgetl(fid);
    quadElements(i,:) = sscanf(tline, '%d')' + 1;
end

% 跳到 # number of geometric entity indices
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, '# number of geometric entity indices')
        break;
    end
end

% 读取数量
numGroupIdx = sscanf(tline, '%d');


while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, '# Geometric entity indices')
        break;
    end
end
groupIdx = zeros(numGroupIdx,1);
for i = 1:numGroupIdx
    tline = fgetl(fid);
    groupIdx(i) = sscanf(tline, '%d');
end

% 分组索引
groups = unique(groupIdx);
index = cell(length(groups),1);
for k = 1:length(groups)
    index{k} = find(groupIdx == groups(k));
end

fclose(fid);
end