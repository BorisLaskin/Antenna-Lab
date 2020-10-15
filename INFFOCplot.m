function STR=INFFOCplot(STR,FullTime)
COLOR=[1 0 0; 0 0.5 0; 0.749 0 0.749;...
    0 0.749 0.749; 0.847 0.161 0; 1 0 1; 0.682 0.467 0; 0.871 0.49 0];
n=numel(STR);
t=zeros(size(STR(1).XZ(:,1)));
for i=1:n
    t=t+STR(i).Time;
end
maxt=max(max(t));

XZnew=zeros(size(STR(n).XZ));
XZold=STR(n).XZ;
Qold=STR(n).Q;
STRnew=STR(n);

dt=FullTime+maxt-t;
XZnew(:,1)=XZold(:,1)+Qold(:,1).*dt(:,1);
XZnew(:,2)=XZold(:,2)+Qold(:,2).*dt(:,1);
STRnew.XZ=XZnew;
STRnew.Time=dt;
STR=[STR STRnew];

hold on;
figure;
if n>0
    NanMass=zeros(size(STR(1).XZ(:,1)));NanMass(:,:)=NaN;
    for i=1:n
    IND=isnan(STR(i+1).XZ(:,2));
    MassDataX=[STR(i).XZ(:,1)';STR(i+1).XZ(:,1)';NanMass'];
    MassDataZ=[STR(i).XZ(:,2)';STR(i+1).XZ(:,2)';NanMass'];
    MassDataX(:,IND)=[];MassDataZ(:,IND)=[];
    DatX=MassDataX(:);
    DatZ=MassDataZ(:);

    h(i,:)=line(DatZ,DatX);
    hold on;
    set(h(i,:),'LineWidth',1,'Color',COLOR(mod(i,8),:));
    end
end
end