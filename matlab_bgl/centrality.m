addpath /home/whuan/.matlab/R2008a/matlab_bgl; %Add directory to search path
clc; %Clear command window
clear all; %Remove all variables, globals, functions and MEX links
t0 = clock; %Current date and time as date vector
mkdir Data; %创建文件夹，储存相关数据文件

%计算数据集里的蛋白质数目(即蛋白质网络中的节点数)
DATAfile = 'Scere20120818_ppi_remove.txt';
fid = fopen(DATAfile,'r'); %只读方式打开
Df = textscan(fid,'%s %s'); %读入文件
fclose(fid); %若关闭成功，返回0，否则返回-1
Pro = union(Df{1},Df{2}); %用union取并集，Pro存放网络中的所有蛋白质
num = length(Pro); %num是蛋白质总数
save Data/Pro.mat Pro; %设置保存的路径，注意，linux下的路径名以斜杠(/)分隔，windows下以反斜杠(\)分隔，但是这里用斜杠(/)在windows下也可以正常运行

%创建邻接矩阵
[bool,edge(:,1)] = ismember(Df{1},Pro); %edge以蛋白质编号的形式存放网络中的所有边
[bool,edge(:,2)] = ismember(Df{2},Pro);
AdjMatrix = zeros(num); %AdjMatrix为邻接矩阵
for i = 1:length(Df{1}) %length(Df{1})是网络的边数
    AdjMatrix(edge(i,1),edge(i,2)) = 1;
    AdjMatrix(edge(i,2),edge(i,1)) = 1;
end
SparseAdjMatrix = sparse(AdjMatrix); %sparse求稀疏矩阵，用于保存，节省空间
save Data/AdjMatrix.mat SparseAdjMatrix;

%计算节点度dc
degree = sum(AdjMatrix,2); %对AdjMatrix的每行求和，degree存放每个节点的度

%计算ic
B = diag(degree)+ones(num)-AdjMatrix; %diag(degree)是以degree为对角线元素，其他元素为零的矩阵
M = pinv(B); %pinv求广义逆矩阵
ic = zeros(num,1); %ic存放每个节点的ic值
for i = 1:num
    for j = 1:num
        if j ~= i
            ic(i) = (M(i,i)+M(j,j)-M(i,j)) + ic(i); %计算ic值
        end
    end
end
ic = num./ic;

%计算ec
[V,D] = eig(AdjMatrix); %对邻接矩阵求其特征值和特征向量；对角矩阵D的对角线上元素是特征值，并且是由小到大排列的，矩阵V的每一列是对应于特征值的特征向量
if sum(V(:,num)<0) < sum(V(:,num)>0)
    ec = V(:,num); %最大的特征值对应的特征向量，即为ec；由于特征值是由小到大排列的，所以最大的特征值在最后一列，即第num列
else
    ec = -V(:,num); %如果特征向量的num个分量中绝大部分都是负值(负值个数大于正值个数)，则对特征向量取负值(即乘以-1)仍然是对应于相应特征值的特征向量，若不取负变正的话，结果很糟糕
end

%计算sc
sc = zeros(num,1); %sc存放每个节点的sc值
for i = 1:num
    for j = 1:num
        sc(i) = exp(D(j,j))*(V(i,j))^2 + sc(i); %计算sc值
    end
end

%计算bc
[bc,E] = betweenness_centrality(SparseAdjMatrix); %matlabBGL中的betweenness_centrality函数；bc返回点介数，E返回边介数

%计算cc
Dist = graphallshortestpaths(SparseAdjMatrix); %graphallshortestpaths函数
Dist(Dist==Inf) = num;
cc = (num-1)./sum(Dist,2); %这里计算的是相对cc值

%保存中心性测度的相关数据
params = [degree,ic,ec,sc,bc,cc];
value = zeros(num,6); pid = zeros(num,6); pname = cell(num,6); sorting = cell(1,6);
for i = 1:6
    [value(:,i),pid(:,i)] = sort(params(:,i),'descend');  %对每种参数降序排列，value保存排序后的值，pid保存排序后的蛋白质编号
    pname(:,i) = Pro(pid(:,i)); %pname存放排序后的蛋白质
    sorting{i} = [num2cell(pid(:,i)),pname(:,i),num2cell(value(:,i))];
end
save Data/params.mat params;
save Data/sorting.mat sorting;

%计算程序运行时间
elapsed_time = etime(clock,t0); %Elapsed time
save Data/elapsed_time.mat elapsed_time;




function seven_Centralities(locs)
clc; %Clear command window
t0 = clock; %Current date and time as date vector
mkdir locs_2015;


%计算数据集里的蛋白质数目(即蛋白质网络中的节点数)
PPI='_PPIW.txt';
DATAfile =strcat(locs,PPI);
fid = fopen(DATAfile,'r'); %只读方式打开
Df = textscan(fid,'%s%s%f64'); %read file
fclose(fid); %若关闭成功，返回0，否则返回-1
Pro = union(Df{1},Df{2}); %用union取并集，Pro存放网络中的所有蛋白质
num = length(Pro); %num是蛋白质总数
DR3=Df{3};
save locs_2015/Pro.mat Pro; %设置保存的路径，注意，linux下的路径名以斜杠(/)分隔，windows下以反斜杠(\)分隔，但是这里用斜杠(/)在windows下也可以正常运行

%创建邻接矩阵
[bool,edge(:,1)] = ismember(Df{1},Pro); %edge以蛋白质编号的形式存放网络中的所有边
[bool,edge(:,2)] = ismember(Df{2},Pro);
AdjMatrix = zeros(num); %AdjMatrix为邻接矩阵
for i = 1:length(Df{1}) %length(Df{1})是网络的边数
    AdjMatrix(edge(i,1),edge(i,2)) = DR3(i);
    AdjMatrix(edge(i,2),edge(i,1)) = DR3(i);
end
SparseAdjMatrix = sparse(AdjMatrix); %sparse求稀疏矩阵，用于保存，节省空间
save locs_2015/AdjMatrix.mat SparseAdjMatrix;

%计算节点度dc
degree = sum(AdjMatrix,2); %对AdjMatrix的每行求和，degree存放每个节点的度

%计算ic
B = diag(degree)+ones(num)-AdjMatrix; %diag(degree)是以degree为对角线元素，其他元素为零的矩阵
M = pinv(B); %pinv求广义逆矩阵
ic = zeros(num,1); %ic存放每个节点的ic值
for i = 1:num
    for j = 1:num
        if j ~= i
            ic(i) = (M(i,i)+M(j,j)-M(i,j)) + ic(i); %计算ic值
        end
    end
end
ic = num./ic;

%计算ec
[V,D] = eig(AdjMatrix); %对邻接矩阵求其特征值和特征向量；对角矩阵D的对角线上元素是特征值，并且是由小到大排列的，矩阵V的每一列是对应于特征值的特征向量
length(V)
if sum(V(:,num)<0) < sum(V(:,num)>0)
    ec = V(:,num); %最大的特征值对应的特征向量，即为ec；由于特征值是由小到大排列的，所以最大的特征值在最后一列，即第num列
else
    ec = -V(:,num); %如果特征向量的num个分量中绝大部分都是负值(负值个数大于正值个数)，则对特征向量取负值(即乘以-1)仍然是对应于相应特征值的特征向量，若不取负变正的话，结果很糟糕
end

%计算sc
sc = zeros(num,1); %sc存放每个节点的sc值
for i = 1:num
    for j = 1:num
        sc(i) = exp(D(j,j))*(V(i,j))^2 + sc(i); %计算sc值
    end
end

%计算bc
[bc,E] = betweenness_centrality(SparseAdjMatrix); %matlabBGL中的betweenness_centrality函数；bc返回点介数，E返回边介数

%计算cc
Dist = graphallshortestpaths(SparseAdjMatrix); %graphallshortestpaths函数
Dist(Dist==Inf) = num;
cc = (num-1)./sum(Dist,2); %这里计算的是相对cc值


%计算 SOECC
ECC= zeros(num,num);  %ECC为边聚集系数矩阵

for i = 1:num
    for j = 1:num
    neighbornumber=0;
     if ( AdjMatrix(i,j) ~=0)
      for k=1:num
       if(k~=j)&&(k~=i)&&(AdjMatrix(i,k)~=0)&&(AdjMatrix(j,k)~=0)
       neighbornumber=neighbornumber+1;
       end
      end
      if((degree(i)>1) && (degree(j)>1))
      ECC(i,j)=neighbornumber/min(sum(AdjMatrix(i,:))-1 ,sum(AdjMatrix(j,:))-1);
      end
     end 
    end
   end
SOECC = sum(ECC,2); %对ECC的每行求和  


%保存中心性测度的相关数据
degree=num2cell(degree); 
ic=num2cell(ic); 
ec=num2cell(ec); 
sc=num2cell(sc); 
bc=num2cell(bc); 
cc=num2cell(cc); 
SOECC=num2cell(SOECC); 

params = [Pro,degree,ic,ec,sc,bc,cc,SOECC];
OUT=strcat('locs_2015/',locs,'.mat') 
save (OUT,'params')