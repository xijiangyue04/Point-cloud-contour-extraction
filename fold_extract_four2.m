% 根据特征分解得到三个特征向量V(:,1)，V(:,2)，V(:,3)。 V(:,1)与V(:,2)两个向量构成平面
%该平面的法向量为V(:,3)  A=V(1,3),B=V(2,3),C=V(3,3), 
%经过的点为mean_neighbor   D=mean_neighbor*V(:,3)
%依据法向量（A,B,C）及平面上的点(x0，y0，z0)得到的平面方程为：A（x-x0）+B（y-y0）+C（z-z0）=0
%即为Ax+By+Cz-D=0  将包围球内三维点投影到平面上，并依据坐标转换将投影到平面的三维点变为二维点
% 依据包围球范围内的平面点，寻找两个最远点并连城直线，计算所有包围球内二维点到直线的距离及不同测来判断折边
%输入变量radius包围球半径,    PL_threshold：点到直线距离中误差限差，值越小越有可能是平面   DP_DS:直线不同侧点数差别几倍
%rank_dis_threshold 第i个点到直线距离排名第几的阈值一般取10 
%距离分辨率为0.005时，radius=0.05,PL_threshold=0.005,DP_DS=4,rank_dis_threshold=3的时候效果比较好
% radius 一般为距离分辨率的10倍，PL_threshold一般与距离分辨率相似，DP_DS=4,rank_dis_threshold=3
function [fold_pnt] = fold_extract_four2(input,radius,PL_threshold,DP_DS, rank_dis_threshold)    %该方法效果比较好
n=size(input,1);
[points] = sphere_points(input,radius);
fold_pnt=[];
for i=1:n

%以下命令是根据主成分或正交整体最小二乘实现包围球平面拟合及各个法向量的计算
    [parameter] = sphere_PCA(points{i}); % 该平面不是TLS平面拟合的平面，而是与TLS_plane垂直的平面
%以下命令是将包围球内三维点投影到平面上，并将其变为一个平面上的三维点
   [project_plane] = PC_Proj(parameter,points{i});
    
% 以下命令是根据坐标转换将投影到平面的三维点变为二维点
    normal1=parameter(:,1:3);
    normal2=[0 0 1];
    [R]=Rotation_matrix(normal1,normal2);
    neighbor_pnt=project_plane'*R;
    x=neighbor_pnt(:,1);
    y=neighbor_pnt(:,2);
      
 %以下命令是依据所有点之间距离寻找最远的那两个点（这几行命令可能有点啰嗦）
    neighbor_number=size(neighbor_pnt,1);
    [neighbor_idx,dis]=knnsearch(neighbor_pnt,neighbor_pnt,'k', neighbor_number); 
    maxdis_idx=find(dis(:,neighbor_number)==max(dis(:,neighbor_number)));  
    farthest_idx=[maxdis_idx,neighbor_idx(maxdis_idx,neighbor_number)];
    farthest_idx= sort(farthest_idx');
    [farthest_idx,~,~]=unique( farthest_idx','rows');

    line_number=size(farthest_idx,1); %确定有几对最远点
    
    %以下命令实现点到最远点连线距离

%     for j=1:1
    farthest_pnt=neighbor_pnt(farthest_idx(1,:),1:2); %两个最远点坐标
    B=[farthest_pnt(:,1),ones(2,1)];L=farthest_pnt(:,2);
    line_parameter=pinv(B'*B)*B'*L;%计算两个最远点连接形成的直线
    PL_dis1=abs(line_parameter(1)*neighbor_pnt(:,1)-neighbor_pnt(:,2)+line_parameter(2))/sqrt(line_parameter(1)^2+1); %所有点到这个直线的距离
    PL_dis=round(PL_dis1*1000)/1000;
%     end
    sort_dis=sort(PL_dis,'descend');%对所有包围球内到直线距离进行排序
     PL_std=std(PL_dis);
   
    %计算包围球内点在直线的左侧及右侧个数
    judge_side=line_parameter(1)*x+line_parameter(2)-y;
    right_side_number=length(find(judge_side>0));
    left_side_number=length(find(judge_side<0));
    if (right_side_number>DP_DS*left_side_number | left_side_number>DP_DS*right_side_number) & PL_std>PL_threshold 
        if length(PL_dis)<rank_dis_threshold+1
         [~,fold_line]=max(PL_dis);
        else
        fold_line=find(PL_dis>sort_dis(rank_dis_threshold+1,:));
        end
        fold_pnt=unique([fold_pnt;points{i}(fold_line,:)],'rows');
    end
end

