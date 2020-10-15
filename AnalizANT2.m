function SREZ=AnalizANT2(STR,TOTrace)
%функци€ весового взвешивани€ Ampl должна включать не только аплитуду луча,
%но и возможность учета пол€ризации
%-------определение направлени€ сканировани€ и точки фокусировки-----------
if TOTrace==1
    %--------фронт волны падает на антенну---------------------------------
    n=numel(STR);
    XZ=STR(n).XZ;
    Q=STR(n).Q;
    CABMain=STR(n).CABMain;
    t=zeros(size(STR(1).XZ(:,1)));
    for i=1:n
        t=t+STR(i).Time;
    end
    maxt=max(max(t));%!!!!!
    
    NaNind=find(isnan(t));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    
    NaNind=find(isnan(XZ(:,2)));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    
    NaNind=find(isnan(Q(:,2)));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    
    dt=maxt-t;
    XZ(:,1)=XZ(:,1)+Q(:,1).*dt(:,1);
    XZ(:,2)=XZ(:,2)+Q(:,2).*dt(:,1);
    
    [dt,~]=fminsearch(@(t)findTime(XZ,Q,t),0);

    
    XZET0=zeros(size(XZ));
    N=numel(XZ(:,2));
    XZET0(:,1)=XZ(:,1)+Q(:,1).*dt;%первое приближение
    XZET0(:,2)=XZ(:,2)+Q(:,2).*dt;%фронта падающей волны
                                    %после последнего отражени€
    T=[sum(XZET0(:,1))/N sum(XZET0(:,2))/N];

    clear('tol','tol1','tol2','NSch','Rold','dT','Told','XYZ','XYZET0','Q');
    clear('flag','flgood','timer','timer2','Time','Mainind','Mng','COSMAIN');
else
    %------------излучение из точки----------------------------------------
    n=numel(STR);
    XZ=STR(n).XZ;
    Q=STR(n).Q;
    CABMain=STR(n).CABMain;
    Ampl=STR(n).Ampl;
    t=zeros(size(STR(1).XZ(:,1)));
    for i=1:n
        t=t+STR(i).Time;
    end
    NaNind=find(isnan(t));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    Ampl(NaNind,:)=[];
    
    NaNind=find(isnan(XZ(:,2)));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    Ampl(NaNind,:)=[];
    
    NaNind=find(isnan(Q(:,2)));
    t(NaNind)=[];
    XZ(NaNind,:)=[];
    Q(NaNind,:)=[];
    CABMain(NaNind,:)=[];
    Ampl(NaNind,:)=[];
    [~,Mainind]=min(abs(CABMain-1));
    Mainind=Mainind(1);
    
    maxt=max(max(t));
    XZET0=zeros(size(XZ));
    XZET0(:,1)=XZ(:,1)+Q(:,1).*(maxt-t(:));%первое приближение
    XZET0(:,2)=XZ(:,2)+Q(:,2).*(maxt-t(:));%фронта падающей волны
                                        %после последнего отражени€
    
    %----------формируем набор отклонений----------------------------------
    proc=0.1; N=3; Alfa=acos(Q(Mainind,1)); Beta=acos(Q(Mainind,2));
    P0=XZET0(Mainind,1)*Q(Mainind,1)+XZET0(Mainind,2)*Q(Mainind,2);
    P=linspace(P0*(1-proc),P0*(1+proc),N);
    Alfa=linspace(Alfa*(1-proc),Alfa*(1+proc),N);
    [MP,MA]=meshgrid(P,Alfa);
    MP=MP(:); MA=MA(:);
    N=numel(MP);
    DSKO=zeros(size(MP)); OPTIM=zeros(N,2);
    fun=@(x)SubfindDSKO(XZET0,x(1),x(2),Ampl);
    options=optimset('tolX',10^-4);
    for i=1:N
        [OPTIM(i,:),DSKO(i)]=fminsearch(fun,[cos(MP(i)) cos(MA(i))],options);
    end
    [~,Mainind]=min(abs(DSKO));
    BestOPTIM=real(OPTIM(Mainind,:));
    %BestM=MSKO(Mainind);
    
    clear('Alfa','Beta','proc','N','MP','MA','MB','fun','Alfa','Beta','P','CABMain','Ampl');
    clear('XYZ','XYZET0','Q','Mainind','COSMAIN','NaNind','P0','t','maxt','MSKO','OPTIM');
end
%--------------------------------------------------------------------------
if TOTrace==1
    n=numel(STR);
    XZ=STR(n).XZ;
    Q=STR(n).Q;
    
    t=zeros(size(STR(1).XZ(:,1)));
    for i=1:n
        t=t+STR(i).Time;
    end
    NaNind=find(isnan(t));
    Ampl=STR(n).Ampl;
    
    BestMOne=SubfindBarC(T,XZ,Q,Ampl,NaNind);
    BarCOne=SubBarC(T,XZ,Q,Ampl);
    BarCOne(NaNind)=NaN;
    fun=@(x)SubfindBarC([x(1) x(2)],XZ,Q,Ampl,NaNind);
    
    options=optimset('tolX',10^-6);
    [OPTIM,BestMTwo]=fminsearch(fun,T,options);
    BarCTwo=SubBarC(OPTIM,XZ,Q,Ampl);
    BarCTwo(NaNind)=NaN;
    
    SREZ=struct('MassOne',BestMOne,'MassDataOne',BarCOne,'FpointOne',T,'MassTwo',BestMTwo,'MassDataTwo',BarCTwo,'FpointTwo',OPTIM);
    
    
    
    %SREZ=struct('MassTwo',BestMTwo,'MassDataTwo',BarCTwo,'FpointTwo',OPTIM);
    clear('BestMOne','BestMTwo','BarCOne','BarCTwo','T','OPTIM','NaNind');
    clear('options','XYZ','Q','Ampl','t');
else
    n=numel(STR);
    XZ=STR(n).XZ;
    Q=STR(n).Q;
    Ampl=STR(n).Ampl;
    
    t=zeros(size(STR(1).XZ(:,1)));
    for i=1:n
        t=t+STR(i).Time;
    end
    NaNind=find(isnan(t));
    
    maxt=max(max(t));
    XZET0=zeros(size(XZ));
    XZET0(:,1)=XZ(:,1)+Q(:,1).*(maxt-t(:));
    XZET0(:,2)=XZ(:,2)+Q(:,2).*(maxt-t(:));
    
    [MSKO,KO]=SubMSKO(XZ,Q,BestOPTIM(1),BestOPTIM(2),Ampl,t);
    %DSKO=SubDSKO(XYZET0,BestOPTIM(1),BestOPTIM(2),BestOPTIM(3),Ampl);
    %DSKO(NaNind)=NaN;
    
    [Zmax,~]=max(XZ(:,2));
    D=-Zmax; 
    dt=zeros(size(STR(1).XZ(:,1)));
    dt(:)=-(XZ(:,2)+D)./Q(:,2);
    XZAPER(:,1)=XZ(:,1)+Q(:,1).*dt(:);
    XZAPER(:,2)=XZ(:,2)+Q(:,2).*dt(:);
    t=t+dt;
    
    SREZ=struct('XZET0',XZET0,'MassData',KO,'Mass',MSKO,'pABC',[BestOPTIM real(sqrt(1-BestOPTIM(2).^2))],'XZAPER',XZAPER,'APERTime',t,'Ampl',Ampl);
    clear('Zmax','Zind','maxt','DSKO','BestM','BestOPTIM','Q','XYZ','D','dt');
    clear('XYZAPER','XYZET0','t');
end
end
function Norm=findTime(XZold,Qold,dt)
    XZnew(:,1)=XZold(:,1)+Qold(:,1).*dt(:,1);
    XZnew(:,2)=XZold(:,2)+Qold(:,2).*dt(:,1);
    ind1=find(isnan(XZnew(:,2)));
    N=numel(XZnew(:,2));
    T=[sum(XZnew(:,1))/N sum(XZnew(:,2))/N];
    
    Norm=sqrt(sum((XZnew(:,1)-T(1)).^2+(XZnew(:,2)-T(2)).^2)/(N-1));
end

function M=SubfindDSKO(XZ,p,A,Ampl) %поискова€ функци€
    DSKO=Ampl(:).^2.*(XZ(:,1)*A+XZ(:,2)*sqrt(1-A^2)-p).^2;
    M=sum(DSKO);
end
function [MSKO,KO]=SubMSKO(XZ,Q,p,A,Ampl,t)
    B=sqrt(1-A^2);
    time=-(A*XZ(:,1)+B*XZ(:,2)-p)./(Q(:,1)*A+Q(:,2)*B);
    FinTime=t+time;
    NaNind=find(isnan(FinTime));
    tTime=FinTime; tTime(NaNind)=[];
    Laverage=sum(tTime)/numel(tTime);
    %KO=(FinTime-Laverage).^2;
    KO=abs(FinTime-Laverage);
    MSKO=sqrt(sum((tTime-Laverage).^2)/(numel(tTime)-1));
    %MSKO=Ampl(:).^2.*(XYZ(:,1)*cos(A)+XYZ(:,2)*cos(B)+XYZ(:,3)*sqrt(1-cos(A)^2-cos(B)^2)-p).^2;
end
function M=SubfindBarC(T,XZ,Q,Ampl,NaNind)%поискова€ функци€
    N=numel(XZ(:,1));
    BarC=zeros(N,1);
    BarC(:)=(T(1)-XZ(:,1)).^2+(T(2)-XZ(:,2)).^2-((T(1)-XZ(:,1)).*Q(:,1)+(T(2)-XZ(:,2)).*Q(:,2)).^2;  
    BarC(NaNind)=[];
    Ampl(NaNind)=[];
    M=sum(Ampl(:).^2.*BarC(:))/sum(Ampl(:).^2);
end
function BarC=SubBarC(T,XZ,Q,Ampl)
    N=numel(XZ(:,1));
    BarC=zeros(N,1);
    BarC(:)=(T(1)-XZ(:,1)).^2+(T(2)-XZ(:,2)).^2-((T(1)-XZ(:,1)).*Q(:,1)+(T(2)-XZ(:,2)).*Q(:,2)).^2;     
end