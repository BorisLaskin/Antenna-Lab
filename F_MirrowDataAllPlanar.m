function  V_Mirrow=F_MirrowDataAllPlanar(V_allParam,SysVal,flagtype)
    %1-st mirrow------------------------------------------
    D=V_allParam.D;
    C=V_allParam.C;
    d=V_allParam.d;
    FOC=V_allParam.FOC;
    delta=V_allParam.delta;
    rad=D/2+C;
%     N1=100;
    N2=1000;
    Q=abs(V_allParam.Q2);
    B=V_allParam.B;
    if flagtype==1

        sgn=sign(B);
    else
        sgn=V_allParam.dfmsgn;

    end
    %---------------------------------------------------------------------
    afamax=asin(rad/FOC);
    afa=linspace(-afamax,afamax,N2);
    f1=FOC;
    f=delta;
    t1=(1-sgn*(f1/d)*sin(afa/2).^2).^(f1/(f1-sgn*d));
    t2=(cos(afa/2).^2).^(d/(d-sgn*f1));
    ro=((1/d)*sin(afa/2).^2+(1/(d-f))*(t1).*t2).^-1;

    Z2=f+ro.*cos(afa);
    X2=ro.*sin(afa);
    Z1=f+ro.*cos(afa)-(d-ro.*sin(afa/2).^2).*(1-(((ro-sgn*f1).*sin(afa))./(2*(d-ro.*sin(afa/2).^2))).^2);
    X1=sgn*f1*sin(afa);
    %---------------------------------------------------------------------
    Tafa=afa;
    TZ2=Z2;
    TX2=X2;
    TZ1=Z1;
    TX1=X1;
    %---------------------------------------------------------------------
    IND=find(X1>D/2+C|X1<-D/2+C);
    afa(IND)=[];
    X1(IND)=[];
    X2(IND)=[];
    Z1(IND)=[];
    Z2(IND)=[];
    AFA=afa;
    %---------------------------------------------------------------------
    [~,IND2] = min(abs(X1-C));
    XC1=X1(IND2);ZC1=Z1(IND2);XC2=X2(IND2);ZC2=Z2(IND2);AFAC=afa(IND2);
    %---------------------------------------------------------------------
    dz1 = sgn*gradient(X1);
    dx1 = -sgn*gradient(Z1);
    a=sqrt(dz1.^2+dx1.^2);
    dz1=dz1./a;dx1=dx1./a;
    %quiver(Z1,X1,dz1,dx1);
    
    %hold on;axis equal;
    %plot(Z1,X1,'r*',Z2,X2,'b*');
    
    dz2 = -gradient(X2);
    dx2 = gradient(Z2);
    a=sqrt(dz2.^2+dx2.^2);
    dz2=dz2./a;dx2=dx2./a;
    %quiver(Z2,X2,dz2,dx2);
    %---------limits-----------------------
    Zmin1=min(Z1);
    Zmin1=min(Zmin1,delta);
    Zmax1=max(Z2);
    Zlength=Zmax1-Zmin1;
    a=381;
    b=424;
    Xmin1=-a*Zlength/b/2+C;
    Xmax1=a*Zlength/b/2+C;
    
    Xmin21=min(X1);Xmin22=min(X2);Xmin2=min([Xmin21,Xmin22]);
    Xmax21=max(X1);Xmax22=max(X2);Xmax2=max([Xmax21,Xmax22]);
    Xlength=Xmax2-Xmin2;
    Zmax2=(Xlength*b+a*Zmin1)/a;
    if Xmin2<Xmin1||Xmax2>Xmax1
        XLIM=[Xmin2 Xmax2];
        ZLIM=[Zmin1 Zmax2];
    else
        XLIM=[Xmin1 Xmax1];
        ZLIM=[Zmin1 Zmax1];
    end
    %---------limits-----------------------
    V_Mirrow(1).X=X1;
    V_Mirrow(1).Z=Z1;
    V_Mirrow(2).X=X2;
    V_Mirrow(2).Z=Z2;
    V_Mirrow(1).dx=dx1;
    V_Mirrow(1).dz=dz1;
    V_Mirrow(2).dx=dx2;
    V_Mirrow(2).dz=dz2;
    V_Mirrow(1).AFA=AFA;
    V_Mirrow(2).AFA=AFA;
    V_Mirrow(1).R=V_allParam.R1;
    V_Mirrow(2).R=V_allParam.R2;
    
    V_Mirrow(1).TX=TX1;
    V_Mirrow(1).TZ=TZ1;
    V_Mirrow(2).TX=TX2;
    V_Mirrow(2).TZ=TZ2;
    V_Mirrow(1).Talfa=Tafa;
    V_Mirrow(2).Talfa=Tafa;
    
    V_Mirrow(1).CX=XC1;
    V_Mirrow(1).CZ=ZC1;
    V_Mirrow(2).CX=XC2;
    V_Mirrow(2).CZ=ZC2;
    V_Mirrow(1).CAFA=AFAC;
    V_Mirrow(2).CAFA=AFAC;
    
    V_Mirrow(1).XLIM=XLIM;
    V_Mirrow(1).ZLIM=ZLIM;
end