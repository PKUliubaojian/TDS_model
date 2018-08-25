function SpecularPoint=SP(PT,PR)
    %PT��PR����N*2�ľ��󣬴���N���㣬���һ��N*2�ľ��󣬴����Ӧ�ľ����
    %��ʼ����������
    SpecularPoint=zeros(size(PT));
    options.Algorithm='interior-point';
    options.Display='off';
    options.OptimalityTolerance=1e-6;
    options.SpecifyObjectiveGradient=true;
    options.UseParallel=false;
    options.HessianFcn='objective';
    sizes=size(PT,1);
    %����ÿ����Եľ����λ��
    for i=1:sizes
        %����ÿ�����䷴����
        Pi=PT(i,:);
        Pr=PR(i,:);
        %�½�Ŀ�꺯�����ݶ�
        newfun=@(f,g,H)SPdis(Pi,Pr,f);
        %��ʼֵ��ѡ��Ϊ���ջ��͵������������򽻵�
        coe=6378137./sqrt(Pr(1).^2+Pr(2).^2+Pr(3).^2/(1-0.00669437999013));
        xstart=coe.*Pr;
        SpecularPoint(i,:)=fmincon(newfun,xstart,[],[],[],[],[],[],'WGS1984_cons',options);
        %[SpecularPoint(i,:),~,~,~]=fmincon(newfun,xstart,[],[],[],[],[],[],'WGS1984_cons',options);
    end
end
