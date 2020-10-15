function STR=INFFOCplotXZ(STR,V_allParam,AXS,Flag,tip)
COLOR=[1 0 0; 0 0.5 0; 0.749 0 0.749;...
    0 0.749 0.749; 0.847 0.161 0; 1 0 1; 0.682 0.467 0; 0.871 0.49 0];
n=numel(STR);
t=zeros(size(STR(1).XYZ(:,1)));
for i=1:n
    t=t+STR(i).Time;
end
maxt=max(max(t));

if tip==1||tip==2
    FullTime=2*abs(V_allParam.S2_delta);
else
    dZ=4*abs(max(STR(2).XYZ(:,3))-min(STR(2).XYZ(:,3)));
    MeshFI=STR(1).MeshFI;
    ind1=find(MeshFI==0);
    ind2=find(MeshFI==pi);
    ind3=find(MeshFI==pi/2);
    ind4=find(MeshFI==3*pi/2);
    ind1=union(ind1,ind2);
    ind2=union(ind3,ind4);
    ind=union(ind1,ind2);
    if tip==4
    FullTime=dZ;
    else
        FullTime=2*abs(V_allParam.S2_delta);
    end
    for i=1:n
        STR(i).XYZ=STR(i).XYZ(ind,:);
        STR(i).Q=STR(i).Q(ind,:);
    end
    t=t(ind,:);
end
XYZnew=zeros(size(STR(n).XYZ));
XYZold=STR(n).XYZ;
Qold=STR(n).Q;
STRnew=STR(n);


dt=FullTime+maxt-t;
XYZnew(:,1)=XYZold(:,1)+Qold(:,1).*dt(:,1);
XYZnew(:,2)=XYZold(:,2)+Qold(:,2).*dt(:,1);
XYZnew(:,3)=XYZold(:,3)+Qold(:,3).*dt(:,1);
STRnew.XYZ=XYZnew;
STRnew.Time=dt;
STR=[STR STRnew];


%figure;
if n>0
    NanMass=zeros(size(STR(1).XYZ(:,1)));NanMass(:,:)=NaN;
    for i=1:n
    
        MassDataX=[STR(i).XYZ(:,1)';STR(i+1).XYZ(:,1)';NanMass'];
        MassDataY=[STR(i).XYZ(:,2)';STR(i+1).XYZ(:,2)';NanMass'];
        MassDataZ=[STR(i).XYZ(:,3)';STR(i+1).XYZ(:,3)';NanMass'];
        IND=isnan(STR(i+1).XYZ(:,3));
        MassDataX(:,IND)=[];MassDataY(:,IND)=[];MassDataZ(:,IND)=[];
        DatX=MassDataX(:);
        DatY=MassDataY(:);
        DatZ=MassDataZ(:);
        if numel(find(IND))<numel(STR(i+1).XYZ(:,3))
            if Flag(i)
                h(i,:)=line(AXS,DatZ,DatX);
                set(h(i,:),'LineWidth',1,'Color',COLOR(mod(i,8),:));
            end
        end
    end
end
end