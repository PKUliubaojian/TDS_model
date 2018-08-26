function [f,g,H]=SPdis(Pt,Pr,PS)
%%�޸�ʱ��2018.8.19��
%����������
%�������߶��ܾ��룺PtPS+PrPS
%%�������
%Pt�����������
%Pr�����ջ�����
%PS��;���ĵ�
%�������
%f��Ϊ���߾���֮��
%g����SP��ƫ����
%h����SP��Hessian����
    f= pdist2(Pt,PS)+pdist2(Pr,PS);
    if nargout>1
        %�����ݶȺ�Hessian����
        g=[4.*PS(1)-2.*Pr(1)-2.*Pt(1);4.*PS(2)-2.*Pr(2)-2.*Pt(2);4.*PS(3)-2.*Pr(3)-2.*Pt(3)];
        H=[4,0,0;0,4,0;0,0,4];
    end
end