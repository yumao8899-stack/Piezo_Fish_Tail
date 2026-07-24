function [voo,zz]= ttq(disp_node , v, nodes_ext)
node_total = size(disp_node,1);
nm = size(v,2);
v_full = zeros(node_total*5, nm);
map=disp_node';
map = map(:); % 拉平成一列
valid = map > 0;
v_full(valid, :) = v(map(valid), :);

voo = nodes_ext;
xvec = v_full;

% 提取各方向位移（假设自由度顺序为x,y,z）
xx = xvec(1:5:end); 
yy = xvec(2:5:end); 
zz = xvec(3:5:end);  

voo(1:node_total , 1) = voo(1:node_total , 1) + xx;
voo(1:node_total , 2) = voo(1:node_total , 2) + yy;
voo(1:node_total , 3) = voo(1:node_total , 3) + zz;