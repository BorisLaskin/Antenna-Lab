function STRout=SimplReflectionEND(STRin,Mir)
%STR-structure content
%STR.XYZ(i,j) i=1,n- control points, j=1,3 -1X;2Y;3Z
%STR.Q(i,j) i=1,n - directive vector, j=1,3 -1X;2Y;3Z
%STR.P1(i,j) i=1,n - polar vector 1, j=1,3 -1X;2Y;3Z
%STR.P2(i,j) i=1,n - polar vector 2, j=1,3 -1X;2Y;3Z
%P1,P2,Q both orthogonal
%Mir.X(Y,Z)- points;
%Mir.ApXY(ApYZ,ApXZ) - aperture in projections;
%--------------------------------------------------------------------------

dx=Mir.dx;dz=Mir.dz;
tx=Mir.X(:);tz=Mir.Z(:);tAFA=Mir.AFA;

tdx=Mir.dx(:);tdz=Mir.dz(:);

MainIND=STRin.MainIND;
MainIND=MainIND(1);
Ampl=STRin.Ampl;
N=numel(Ampl);

XZold=STRin.XZ;
Qold=STRin.Q;
OldNaNind=find(isnan(XZold(:,2)));

XZnew=zeros(size(STRin.XZ));
NORM=zeros(size(STRin.XZ));
Qnew=zeros(size(STRin.XZ));

CX=Mir.CX;
CZ=Mir.CZ;
%-------------------------------------------------------------------------

[OLDX,ttx]=meshgrid(XZold(:,1),tx);
[OLDZ,ttz]=meshgrid(XZold(:,2),tz);

TQNEWX=ttx-OLDX;
TQNEWZ=ttz-OLDZ;

TQoldX=repmat(Qold(:,1)',[numel(tx) 1]);
TQoldZ=repmat(Qold(:,2)',[numel(tx) 1]);

COSQQ=(TQNEWX.*TQoldX+TQNEWZ.*TQoldZ)./sqrt(TQNEWX.^2+TQNEWZ.^2)./sqrt(TQoldX.^2+TQoldZ.^2);

[~,ind]=max(COSQQ,[],1);
XZnew(:,1)=tx(ind(:));
XZnew(:,2)=tz(ind(:));

NORM(:,1)=tdx(ind(:));
NORM(:,2)=tdz(ind(:));
AFA=tAFA(ind(:));
dalfa=abs(tAFA(1)-tAFA(10));
%-------------------------------------------------------------------------         

tol=10^-9;%точность определения места отражения
maxdelta=2*abs(tol);
clear('t','P1Z','P2Z','ind1','ind2');
%-------------------------------------------------------------------------
TEMP=XZnew;
%-------------------------------------------------------------------------
chet=0;

FX=@(alfa)interp1(tAFA,tx,alfa,'spline');
FZ=@(alfa)interp1(tAFA,tz,alfa,'spline');
FdX=@(alfa)interp1(tAFA,tdx,alfa,'spline');
FdZ=@(alfa)interp1(tAFA,tdz,alfa,'spline');

while abs(maxdelta)>tol&&chet<100
    %------------обнуляем переменные
    MirXTemp=[];MirZTemp=[];MirAFATemp=[];
    
    for i=1:N
        talfa=linspace(AFA(i)-dalfa,AFA(i)+dalfa,100);
        MirXTemp=[MirXTemp FX(talfa)];
        MirZTemp=[MirZTemp FZ(talfa)];
        MirAFATemp=[MirAFATemp talfa];
    end
    OLDX=[];OLDZ=[];ttx=[];ttz=[];TQNEWX=[];TQNEWZ=[];TQoldX=[];TQoldZ=[];COSQQ=[];ind1=[];
    [OLDX,ttx]=meshgrid(XZold(:,1),MirXTemp);
    [OLDZ,ttz]=meshgrid(XZold(:,2),MirZTemp);
    
    TQNEWX=ttx-OLDX;
    TQNEWZ=ttz-OLDZ;

    TQoldX=repmat(Qold(:,1)',[numel(MirXTemp) 1]);
    TQoldZ=repmat(Qold(:,2)',[numel(MirZTemp) 1]);
    
    COSQQ=(TQNEWX.*TQoldX+TQNEWZ.*TQoldZ)./sqrt(TQNEWX.^2+TQNEWZ.^2)./sqrt(TQoldX.^2+TQoldZ.^2);
    [CosVar,ind1]=max(COSQQ,[],1);
    
    XZnew(:,1)=MirXTemp(ind1(:));
    XZnew(:,2)=MirZTemp(ind1(:));
    maxdelta=max(max(abs(TEMP(:,2)-XZnew(:,2))));
    dalfa=dalfa/10;
    chet=chet+1;
end
NORM(:,1)=FdX(MirAFATemp(ind1));
NORM(:,2)=FdZ(MirAFATemp(ind1));
%--------------------------------------------------------------------------

Qnew(:,1)=Qold(:,1)-2*NORM(:,1).*(NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2));
Qnew(:,2)=Qold(:,2)-2*NORM(:,2).*(NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2));

CABMain=zeros(numel(Ampl),1);
CABMain(:)=Qnew(MainIND,1).*Qnew(:,1)+Qnew(MainIND,2).*Qnew(:,2);
time=sqrt((XZnew(:,1)-XZold(:,1)).^2+(XZnew(:,2)-XZold(:,2)).^2);

NaNind=find(abs(1-CosVar)>0.0001);
CABMain(NaNind)=NaN;
time(NaNind)=NaN;
Qnew(NaNind,:)=NaN;
XZnew(NaNind,:)=NaN;
%OldNaNind
CABMain(OldNaNind)=NaN;
time(OldNaNind)=NaN;
Qnew(OldNaNind,:)=NaN;
XZnew(OldNaNind,:)=NaN;
%-----------------------поиск неправильно отраженных от внешней нормали
R=[0 Mir.R];
CMR=Qold(MainIND,1).*R(1,1)+Qold(MainIND,2).*R(1,2);
CNQnew=NORM(:,1).*Qnew(:,1)+NORM(:,2).*Qnew(:,2);
CNQold=NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2);
SignQnewQold=sign(CNQold(:,1))-sign(CNQnew(:,1));
if CMR>=0
    ind=find(SignQnewQold==-2);
else
    ind=find(SignQnewQold==2);
end
Qnew(ind,:)=NaN;
XZnew(ind,:)=NaN;
time(ind)=NaN;
CABMain(ind)=NaN;

STRout=struct('XZ',XZnew,'Q',Qnew,'Time',time,'MainIND',MainIND,'CABMain',CABMain,'Ampl',Ampl);

end

% dx=Mir.dx;dz=Mir.dz;
% tx=Mir.X(:);tz=Mir.Z(:);
% 
% tdx=Mir.dx(:);tdz=Mir.dz(:);
% 
% MainIND=STRin.MainIND;
% MainIND=MainIND(1);
% Ampl=STRin.Ampl;
% 
% XZold=STRin.XZ;
% Qold=STRin.Q;
% 
% XZnew=zeros(size(STRin.XZ));
% NORM=zeros(size(STRin.XZ));
% Qnew=zeros(size(STRin.XZ));
% 
% CX=Mir.CX;
% CZ=Mir.CZ;
% %-------------------------------------------------------------------------
% 
% [OLDX,ttx]=meshgrid(XZold(:,1),tx);
% [OLDZ,ttz]=meshgrid(XZold(:,2),tz);
% 
% TQNEWX=ttx-OLDX;
% TQNEWZ=ttz-OLDZ;
% 
% TQoldX=repmat(Qold(:,1)',[numel(tx) 1]);
% TQoldZ=repmat(Qold(:,2)',[numel(tx) 1]);
% 
% COSQQ=(TQNEWX.*TQoldX+TQNEWZ.*TQoldZ)./sqrt(TQNEWX.^2+TQNEWZ.^2)./sqrt(TQoldX.^2+TQoldZ.^2);
% 
% [~,ind]=max(COSQQ,[],1);
% XZnew(:,1)=tx(ind(:));
% XZnew(:,2)=tz(ind(:));
% 
% NORM(:,1)=tdx(ind(:));
% NORM(:,2)=tdz(ind(:));
% 
% %-------------------------------------------------------------------------         
% 
% tol=10^-9;%точность определения места отражения
% maxdelta=2*abs(tol);
% clear('t','P1Z','P2Z','ind1','ind2');
% %-------------------------------------------------------------------------
% TEMP=XZnew;
% %-------------------------------------------------------------------------
% chet=0;
% while abs(maxdelta)>tol&&chet<100
%     %------------обнуляем переменные
%     A=zeros(size(XZnew(:,1)));
%     B=zeros(size(XZnew(:,1)));
%     C=zeros(size(XZnew(:,1)));
%     time=zeros(size(XZnew(:,1)));
%     
%     %-основной алгоритм
%     indold=isnan(XZold(:,2));
%     if numel(find(indold))>0  
%         XZold(indold,2)=0;%сбоить может тут
%     end
%     
%     A(:,1)=NORM(:,1);
%     B(:,1)=NORM(:,2);
%     C(:,1)=-A(:,1).*XZnew(:,1)-B(:,1).*XZnew(:,2);
%     
%     timeZnam(:,1)=(A(:,1).*Qold(:,1)+B(:,1).*Qold(:,2));
%     time(:,1)=-(A(:,1).*XZold(:,1)+B(:,1).*XZold(:,2)+C(:,1))./timeZnam(:,1);
%     XZnew(:,1)=XZold(:,1)+Qold(:,1).*time(:,1);
%     %----ошибка вычисления пересечения
%     INDNAN=find(timeZnam==0);
%     NFlagNaN=numel(INDNAN);
%     if NFlagNaN
%         XZnew(INDNAN,1)=TEMP(INDNAN,1);
%     end
%     if NFlagNaN
%         XZnew(INDNAN,2)=NaN;
%     end
%     
%     XZnew(:,2)=interp1(tx,tz,XZnew(:,1));
%     maxdelta=max(max(abs(TEMP(:,2)-XZnew(:,2)))); % и тут тоже может сбоить
%     
%     if numel(find(isnan(XZnew(:,2))))==numel(XZnew(:,2))
%         maxdelta=2*tol;
%     end
%     
%     NORM(:,1)=interp1(tx,tdx,XZnew(:,1));
%     NORM(:,2)=interp1(tx,tdz,XZnew(:,1));
%     TEMP=XZnew;
%     
%     NaNXYZind=find(isnan(XZnew(:,2)));
% 
%     chet=chet+1;
% end
% %----------------------------устранение ошибок
% if NaNXYZind
%     XZnew(NaNXYZind,2)=NaN;
%     NORM(NaNXYZind,1)=NaN;
%     NORM(NaNXYZind,2)=NaN;
% end
% 
% 
% clear('INDNAN','NFlagNaN');
% clear('i','A','B','C','D','TEMP','I','J','indold','indnew');
% %--------------------------------------------------------------------------
% 
% Qnew(:,1)=Qold(:,1)-2*NORM(:,1).*(NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2));
% Qnew(:,2)=Qold(:,2)-2*NORM(:,2).*(NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2));
% 
% CABMain=zeros(numel(time),1);
% CABMain(:)=Qnew(MainIND,1).*Qnew(:,1)+Qnew(MainIND,2).*Qnew(:,2);
% 
% CABMain(NaNXYZind)=NaN;
% time(NaNXYZind)=NaN;
% Qnew(NaNXYZind,:)=NaN;
% 
% XZnew(NaNXYZind,:)=NaN;
% %-----------------------поиск неправильно отраженных от внешней нормали
% R=[0 Mir.R];
% CMR=Qold(MainIND,1).*R(1,1)+Qold(MainIND,2).*R(1,2);
% CNQnew=NORM(:,1).*Qnew(:,1)+NORM(:,2).*Qnew(:,2);
% CNQold=NORM(:,1).*Qold(:,1)+NORM(:,2).*Qold(:,2);
% SignQnewQold=sign(CNQold(:,1))-sign(CNQnew(:,1));
% if CMR>=0
%     ind=find(SignQnewQold==-2);
% else
%     ind=find(SignQnewQold==2);
% end
% Qnew(ind,:)=NaN;
% XZnew(ind,:)=NaN;
% time(ind)=NaN;
% CABMain(ind)=NaN;
% 
% STRout=struct('XZ',XZnew,'Q',Qnew,'Time',time,'MainIND',MainIND,'CABMain',CABMain,'Ampl',Ampl);