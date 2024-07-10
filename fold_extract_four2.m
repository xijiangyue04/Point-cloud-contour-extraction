% ���������ֽ�õ�������������V(:,1)��V(:,2)��V(:,3)�� V(:,1)��V(:,2)������������ƽ��
%��ƽ��ķ�����ΪV(:,3)  A=V(1,3),B=V(2,3),C=V(3,3), 
%�����ĵ�Ϊmean_neighbor   D=mean_neighbor*V(:,3)
%���ݷ�������A,B,C����ƽ���ϵĵ�(x0��y0��z0)�õ���ƽ�淽��Ϊ��A��x-x0��+B��y-y0��+C��z-z0��=0
%��ΪAx+By+Cz-D=0  ����Χ������ά��ͶӰ��ƽ���ϣ�����������ת����ͶӰ��ƽ�����ά���Ϊ��ά��
% ���ݰ�Χ��Χ�ڵ�ƽ��㣬Ѱ��������Զ�㲢����ֱ�ߣ��������а�Χ���ڶ�ά�㵽ֱ�ߵľ��뼰��ͬ�����ж��۱�
%�������radius��Χ��뾶,    PL_threshold���㵽ֱ�߾���������޲ֵԽСԽ�п�����ƽ��   DP_DS:ֱ�߲�ͬ�������𼸱�
%rank_dis_threshold ��i���㵽ֱ�߾��������ڼ�����ֵһ��ȡ10 
%����ֱ���Ϊ0.005ʱ��radius=0.05,PL_threshold=0.005,DP_DS=4,rank_dis_threshold=3��ʱ��Ч���ȽϺ�
% radius һ��Ϊ����ֱ��ʵ�10����PL_thresholdһ�������ֱ������ƣ�DP_DS=4,rank_dis_threshold=3
function [fold_pnt] = fold_extract_four2(input,radius,PL_threshold,DP_DS, rank_dis_threshold)    %�÷���Ч���ȽϺ�
n=size(input,1);
[points] = sphere_points(input,radius);
fold_pnt=[];
for i=1:n

%���������Ǹ������ɷֻ�����������С����ʵ�ְ�Χ��ƽ����ϼ������������ļ���
    [parameter] = sphere_PCA(points{i}); % ��ƽ�治��TLSƽ����ϵ�ƽ�棬������TLS_plane��ֱ��ƽ��
%���������ǽ���Χ������ά��ͶӰ��ƽ���ϣ��������Ϊһ��ƽ���ϵ���ά��
   [project_plane] = PC_Proj(parameter,points{i});
    
% ���������Ǹ�������ת����ͶӰ��ƽ�����ά���Ϊ��ά��
    normal1=parameter(:,1:3);
    normal2=[0 0 1];
    [R]=Rotation_matrix(normal1,normal2);
    neighbor_pnt=project_plane'*R;
    x=neighbor_pnt(:,1);
    y=neighbor_pnt(:,2);
      
 %�����������������е�֮�����Ѱ����Զ���������㣨�⼸����������еㆪ�£�
    neighbor_number=size(neighbor_pnt,1);
    [neighbor_idx,dis]=knnsearch(neighbor_pnt,neighbor_pnt,'k', neighbor_number); 
    maxdis_idx=find(dis(:,neighbor_number)==max(dis(:,neighbor_number)));  
    farthest_idx=[maxdis_idx,neighbor_idx(maxdis_idx,neighbor_number)];
    farthest_idx= sort(farthest_idx');
    [farthest_idx,~,~]=unique( farthest_idx','rows');

    line_number=size(farthest_idx,1); %ȷ���м�����Զ��
    
    %��������ʵ�ֵ㵽��Զ�����߾���

%     for j=1:1
    farthest_pnt=neighbor_pnt(farthest_idx(1,:),1:2); %������Զ������
    B=[farthest_pnt(:,1),ones(2,1)];L=farthest_pnt(:,2);
    line_parameter=pinv(B'*B)*B'*L;%����������Զ�������γɵ�ֱ��
    PL_dis1=abs(line_parameter(1)*neighbor_pnt(:,1)-neighbor_pnt(:,2)+line_parameter(2))/sqrt(line_parameter(1)^2+1); %���е㵽���ֱ�ߵľ���
    PL_dis=round(PL_dis1*1000)/1000;
%     end
    sort_dis=sort(PL_dis,'descend');%�����а�Χ���ڵ�ֱ�߾����������
     PL_std=std(PL_dis);
   
    %�����Χ���ڵ���ֱ�ߵ���༰�Ҳ����
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

