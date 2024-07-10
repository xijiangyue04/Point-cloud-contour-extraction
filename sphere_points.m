%�������ݰ�Χ���еĵ������ж��Ƿ�������
%�������input(nx3)  ��İ뾶radius
%�������number_points ���еĵ���
function [points] = sphere_points(input,radius)
n=size(input,1);
for i=1:n
  %�����ó���Ϊradius�����η���ȷ�����·�Χ
line_xyz=find(input(:,1)>input(i,1)-radius & input(:,1)<input(i,1)+radius & input(:,2)>input(i,2)-radius & input(:,2)<input(i,2)+radius & input(:,3)>input(i,3)-radius & input(:,3)<input(i,3)+radius);
input_xyz=input(line_xyz,:);
k=size(input_xyz,1);
 D=sqrt(sum((input_xyz-repmat(input(i,:),k,1)).^2,2));
 sphere_idx=find(D<radius);
 points{i}=input(line_xyz(sphere_idx),:);
end

for i=1:n
    vector{i}=points{i}-input(i,:);
    vector{i}(1,:)=[];
end

    