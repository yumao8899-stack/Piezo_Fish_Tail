function [nodes_ext1,zz]= ttd(disp_node , v, nodes_ext,z_scale_factor)
nodes_ext=nodes_ext*1e3;
nodes_ext(:,1)=nodes_ext(:,1)-110;
node_total = size(disp_node,1);
v_full = zeros(node_total*5, 1);
map=disp_node';
map = map(:); % 拉平成一列
valid = map > 0;
v_full(valid, :) = v(map(valid), :);

% 提取各方向位移（假设自由度顺序为x,y,z）
xx = v_full(1:5:end); 
yy = v_full(2:5:end); 
zz = v_full(3:5:end);  

nodes_ext1(1:node_total , 1) = nodes_ext(1:node_total , 1) + xx;
nodes_ext1(1:node_total , 2) = nodes_ext(1:node_total , 2) + yy;
nodes_ext1(1:node_total , 3) = nodes_ext(1:node_total , 3) + zz* z_scale_factor;
