function [Res,Errorflg] = GetDataReflection(V_Mirrow,V_Angle,V_param)
Errorflg=0;

for i=1:numel(V_Mirrow)
    NewMirrow(i)=V_Mirrow(i); %#ok<AGROW>
    NewMirrow2(numel(V_Mirrow)-i+1)=NewMirrow(i); %#ok<AGROW>
end
Level=0.5;
TOTrace=1; % - тип прокладываемой трассировки 1- В; 2- ИЗ антенны
%------------------Инициализация исходных данных расчетной части-----------
sTH=V_Angle;
sTH=sTH/180*pi;
Nr=25;% число лучей
r=[linspace(-V_param.D/2,V_param.D/2,Nr),0];
r=unique(r);
N=numel(sTH);
ALLDATA1=cell(N,2);
%------------------------------сагиттальная пл-ть--------------------------
h=waitbar(0,'Часть 1 ... ждите...');
clear('STR','STRR','SREZ');
for i=1:N
    S1=struct('r',r,'TH',sTH(i));
    STR=MirReflection2(NewMirrow,S1,S1,TOTrace,Level);
    try
        SREZ=AnalizANT2(STR,TOTrace); %#ok<*AGROW,*SAGROW>
        ALLDATA1{i,1}=STR;
        ALLDATA1{i,2}=SREZ;
    catch
        Errorflg=1;
    end
    waitbar(i/N,h);
end
close(h);
dmass1=zeros(N,2);
try
    for i=1:N
        dmass1(i,:)=ALLDATA1{i,2}.FpointTwo;
    end
catch
    Errorflg=1;
end
dmasstrue=dmass1;
if ~Errorflg
    TOTrace=2;
    ALLDATA2=cell(N,2);
    
    
    h=waitbar(0,'Часть 2 ... ждите...');
    clear('STR','STRR','SREZ');
    for i=1:N
        S2=struct('dx',dmasstrue(i,1),'dz',dmasstrue(i,2));
        STR=MirReflection2(NewMirrow2,S2,S2,TOTrace,Level);
        try
            SREZ=AnalizANT2(STR,TOTrace); %#ok<*AGROW,*SAGROW>
            ALLDATA2{i,1}=STR;
            ALLDATA2{i,2}=SREZ;
        catch
            Errorflg=1;
        end
        waitbar(i/N,h);
    end
    Res=struct('Sagit',ALLDATA1,'OutTrace',ALLDATA2,'Dmass',dmasstrue);
    close(h);
else
    Res=struct('Sagit',ALLDATA1);
end
%
end
