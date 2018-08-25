function SpecularPoint=SP(PT,PR)
    %PT和PR都是N*2的矩阵，代表N个点，输出一个N*2的矩阵，代表对应的镜面点
    %初始化镜面点矩阵
    SpecularPoint=zeros(size(PT));
    options.Algorithm='interior-point';
    options.Display='off';
    options.OptimalityTolerance=1e-6;
    options.SpecifyObjectiveGradient=true;
    options.UseParallel=false;
    options.HessianFcn='objective';
    sizes=size(PT,1);
    %计算每个点对的镜面点位置
    for i=1:sizes
        %迭代每个入射反射点对
        Pi=PT(i,:);
        Pr=PR(i,:);
        %新建目标函数和梯度
        newfun=@(f,g,H)SPdis(Pi,Pr,f);
        %初始值点选择为接收机和地心连线与椭球交点
        coe=6378137./sqrt(Pr(1).^2+Pr(2).^2+Pr(3).^2/(1-0.00669437999013));
        xstart=coe.*Pr;
        SpecularPoint(i,:)=fmincon(newfun,xstart,[],[],[],[],[],[],'WGS1984_cons',options);
        %[SpecularPoint(i,:),~,~,~]=fmincon(newfun,xstart,[],[],[],[],[],[],'WGS1984_cons',options);
    end
end
