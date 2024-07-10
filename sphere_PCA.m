
%采用的是正交总体最小二乘的形式：A(x-xj)+B(y-yj)+C(z-zj)=0  Ax+By+Cz+D=0 （A,B,C）向量为最大特征值所对应的特征向量
%输入变量:input_pnts(nx3)   

function [parameter] = sphere_PCA(input_pnts)
   mean_neighbor=mean(input_pnts,1);
   difference=input_pnts-mean_neighbor;
    M=difference'*difference;
    [V, lamda]=svd(M);
    A=V(1,1);B=V(2,1);C=V(3,1); % 该平面不是TLS平面拟合的平面，而是与TLS_plane垂直的平面
    D=mean_neighbor*V(:,1);
    parameter=[A,B,C,D];