function [Res,Errorflg] = GetDataReflection(V_Mirrow,V_Angle,V_param,AngleMerid)


fi=[0 pi];

numF=1;
ALLDATA1=cell(numF,N,3);
h=waitbar(0,'„асть 1 ... ждите...');
for j=1:numF   
    
    clear('STR','STRR','SREZ');
    for i=1:N
        S1=struct('r',r,'fi',fi,'TH',TH(i),'FI',FI(i));

        %------------------«апуск расчетной части   -------------------------------
        [STR,STRR]=MirReflection2(NewMirrow,S1,S1,TOTrace,Level);
        try
        SREZ(i)=AnalizANT2(STR,TOTrace); %#ok<*AGROW,*SAGROW>
        ALLDATA1{j,i,3}=SREZ(i);
        catch
            subErrorflg(1)=1;
        end
        ALLDATA1{j,i,1}=STR; ALLDATA1{j,i,2}=STRR;
        waitbar(((j-1)*N+i)/numF/N,h);
    end
   
end
clear('STR','STRR');
close(h);

%-------------------------------------------------------------------------
dmass1=zeros(N,3);
dmass2=zeros(N,3);
dmass3=zeros(N,3);
try
    for i=1:N
        dmass1(i,:)=ALLDATA1{j,i,3}.FpointTwo;
    end
catch
    subErrorflg(1)=1;
end
try
    for i=1:N
        dmass2(i,:)=ALLDATA2{j,i,3}.FpointTwo;
    end
catch
    subErrorflg(2)=1;
end
for i=1:N
    dmass3(i,:)=ALLDATA3{j,i,3}.FpointTwo;
end
   
    
%dmasstrue=dmass3;%÷ентр масс кольцевых точек
dmasstrue=dmass3;

if subErrorflg(1)&&subErrorflg(2)&&subErrorflg(3)
    Errorflg=1;
end
if subErrorflg(1)
    dmass1=dmass3;
end
if subErrorflg(2)
    dmass1=dmass3;
end

    
if ~Errorflg
    TOTrace=2;
    ALLDATA4=cell(numF,N,3);
    THL=0;THT=0; %  отклонение оси излучени€ рупора от центра раскрыва 
    % L -  вертикальна€ плоскость ; T - горизонтальна€ поверхность
    th=[0:pi/32:pi/2];
    
    %th=[0:pi/64:pi/2];
    fi=[0:pi/16:2*pi];
    h=waitbar(0,'„асть 4 ... ждите...');
    [~,ZeroInd]=min(abs(V_Angle));
    
    for j=1:numF      

        clear('STR2','STRR2','SREZ2');

        for i=1:N
            S2=struct('dx',dmasstrue(i,1),'dy',dmasstrue(i,2),'dz',dmasstrue(i,3),'THL',THL,'THT',THT,'th',th','fi',fi);

            %------------------«апуск расчетной части   -------------------------------
            [STR2,STRR2]=MirReflection2(NewMirrow2,S2,S2,TOTrace,Level);
            try
                if i==ZeroInd
                    STR2=AnalizANToutside(STR2,V_param.B);
                end
                SREZ2(i)=AnalizANT2(STR2,TOTrace);
                ALLDATA4{j,i,3}=SREZ2(i);
            catch
            end
            ALLDATA4{j,i,1}=STR2; ALLDATA4{j,i,2}=STRR2; 
            waitbar(((j-1)*N+i)/numF/N);
        end
        %------------------------------------------------------------------
        %-----расчет траектории центрального луча
    end
    %---------------------------------------------------------------------
    %S2=struct('dx',dmasstrue(i,1),'dy',dmasstrue(i,2),'dz',dmasstrue(i,3),'THL',THL,'THT',THT,'th',th','fi',fi);
    %---------------------------------------------------------------------
    Res=struct('Sagit',ALLDATA1,'Merid',ALLDATA2,'Full',ALLDATA3,'OutTrace',ALLDATA4,'Dmass',struct('sag',dmass1,'mer',dmass2,'true',dmasstrue));
    close(h);
else
    Res=struct('Sagit',ALLDATA1,'Merid',ALLDATA2,'Full',ALLDATA3);
    Errorflg=0;
end




clear('THL','THT','fi','th','SREZ2','STR2','STRR2','point');

%clear('S1','S2','h','dmass');
    %mesh(Mir.Y,Mir.Z,Mir.X,ones(size(Mir.X)));
save('Antdata','-mat');
clear('ALLDATA','ALLDATA1','ALLDATA2','F','FI','JX','Mir','N','R0','S1','S2','TH','TOTrace');
clear('dmass','dmass1','dmasstrue','dx0','dy0','h','i','j','numF','str');
end
