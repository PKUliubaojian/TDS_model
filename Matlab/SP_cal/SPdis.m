function [f,g,H]=SPdis(Pt,Pr,PS)
%%修改时间2018.8.19日
%作者刘宝剑
%计算折线段总距离：PtPS+PrPS
%%输入参数
%Pt：发射机坐标
%Pr：接收机坐标
%PS：途径的点
%输出参数
%f：为折线距离之和
%g：对SP求偏导数
%h：对SP求Hessian矩阵
    f= pdist2(Pt,PS)+pdist2(Pr,PS);
    if nargout>1
        %计算梯度和Hessian矩阵
        %这里的jac和Hess计算都有误,2019/3/3
        g=[4.*PS(1)-2.*Pr(1)-2.*Pt(1);4.*PS(2)-2.*Pr(2)-2.*Pt(2);4.*PS(3)-2.*Pr(3)-2.*Pt(3)];
        H=[4,0,0;0,4,0;0,0,4];
    end
end