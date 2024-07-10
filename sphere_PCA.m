
%���õ�������������С���˵���ʽ��A(x-xj)+B(y-yj)+C(z-zj)=0  Ax+By+Cz+D=0 ��A,B,C������Ϊ�������ֵ����Ӧ����������
%�������:input_pnts(nx3)   

function [parameter] = sphere_PCA(input_pnts)
   mean_neighbor=mean(input_pnts,1);
   difference=input_pnts-mean_neighbor;
    M=difference'*difference;
    [V, lamda]=svd(M);
    A=V(1,1);B=V(2,1);C=V(3,1); % ��ƽ�治��TLSƽ����ϵ�ƽ�棬������TLS_plane��ֱ��ƽ��
    D=mean_neighbor*V(:,1);
    parameter=[A,B,C,D];