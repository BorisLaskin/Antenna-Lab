function STR=MirReflection2(Mir,S1,S2,TOTrace,Level)

    n=numel(Mir);
    if n>0
        switch TOTrace
%----------------------------Случай 1--------------------------------------
            case 1
            r=S1.r;
            TH=S1.TH;
            %------------------------------Ishodnij front------------------
            CX=Mir(1).CX;%центр апертуры
            CZ=Mir(1).CZ;
            %ДОДЕЛАТЬ!!!!
            L0=2*abs(max(max(Mir(1).Z))-min(min(Mir(1).Z)));% - определение минимальной величины отклонения фронта %ДОДЕЛАТЬ!!!!
             
            %-----------------определение центра апертуры входного зрачка--
            A=-sin(TH);C=-cos(TH); %- вектор нормали
            X0=CX-A*L0;
            Z0=CZ-C*L0;
            %-----------------формирование входного пучка------------------

            Q=zeros(numel(r),2);
            Q(:,1)=A; Q(:,2)=C; Q=Q/sqrt(A^2+C^2);

            %----------------------------------------------end-------------
            XZ=zeros(numel(r),2);
            XZ(:,1)=X0+r*cos(TH);
            XZ(:,2)=Z0-r*sin(TH);

            %-----------------Задание вспомогательных полей----------------
            Nanaliz=numel(XZ(:,1));
            t=zeros(Nanaliz,1);
            
            %Ampl=ones(Nanaliz,1);
            %AmplR=ones(Nanaliz,1);
            k=acos(Level)/max(r);
            Ampl=cos(r*k);
            [~,MainIND]=min(abs(r));
            
            CABMain=ones(numel(r),1);%kosinus ugla megdu gklavnim luchum 
            %------------------------------------------------------
            
            STRin=struct('XZ',XZ,'Q',Q,'Time',t,'MainIND',MainIND,'CABMain',CABMain,'Ampl',Ampl);
            
            clear('A','C');
            clear('r','fi','tmeshR','tmeshFI','L0','TH','FI','X0','Y0','Z0','Nanaliz');
            clear('MeshFI','MeshR','XYZR','MeshFIR','MeshRR','MainINDR');
            
%----------------------------Случай 2--------------------------------------            
            case 2
            %-----------------формирование входного пучка------------------
            %------------------------центральный луч-----------------------
            DX=S2.dx;
            DZ=S2.dz;
            
            CX=Mir(1).CX;%центр апертуры
            CZ=Mir(1).CZ;
            Q0=[(CX-DX),(CZ-DZ)];
            PSIB=atan(Q0(1)/Q0(2));
           
            %---формирование лучей в СК главного луча----------------------
            
            Q=zeros(numel(Mir(1).X),2);
            Q(:,1)=Mir(1).X(:)-DX; Q(:,2)=Mir(1).Z(:)-DZ; 
            Q=Q./sqrt(Q(:,1).^2+Q(:,2).^2);
            thRupor=acos(Q(:,1).*sin(PSIB)+Q(:,2).*cos(PSIB));
            ThMax=max(thRupor);% определение максимального отклонения луча
            %----------------------------------------            
            
            N=25;th=linspace(-ThMax,ThMax,N)+PSIB; th=th(:);
            clear('Q','thRupor');
            Q=zeros(N,2);Q(:,1)=sin(th);Q(:,2)=cos(th);Q=Q./sqrt(Q(:,1).^2+Q(:,2).^2);
            thRupor=acos(Q(:,1).*sin(PSIB)+Q(:,2).*cos(PSIB));
            
            XZ=zeros(N,2);
            XZ(:,1)=XZ(:,1)+DX;
            XZ(:,2)=XZ(:,2)+DZ;

            %-------------------------- end -------------------------------
            %-----------------Задание вспомогательных полей----------------
            Nanaliz=numel(XZ(:,1));
            t=zeros(Nanaliz,1);
            Ampl=RuporTestAmpl(Level,thRupor,ThMax);
            
            
            [~,MainIND]=min(abs(thRupor)); %потенциальная ошибка
                        
            CABMain=ones(Nanaliz,1);
            CABMain(:)=cos(thRupor);%kosinus ugla megdu gklavnim luchum 
            %------------------------------------------------------
            STRin=struct('XZ',XZ,'Q',Q,'Time',t,'MainIND',MainIND,'CABMain',CABMain,'Ampl',Ampl);
               
            clear('n1','MeshFI','MeshTH','m','A','B','C','IA','Nanaliz');
            otherwise
        end
        clear('XYZ','XYZR','Q','QR','P1','P2','Time','MainIND','CABMain','Ampl');
        
        STR(1)=STRin;
        %---------block obrabotki kagdoj poverhnosti
        for i=1:n
            STRout=SimplReflectionEND(STRin,Mir(i));
            STR(i+1)=STRout; %#ok<*AGROW>
            STRin=STRout;
        end
        clear('STRin','STRout');
        if TOTrace==1
            clear('STRinR','STRoutR');
        end
    end

end
function res=RuporTestAmpl(Level,th,ThMax)
    J0=@(x)1/pi*integral(@(t)cos(x.*sin(t)),0,pi);
    a=linspace(0,100,1000);
    f1=zeros(size(a));
    for i=1:numel(a)
        f1(i)=J0(a(i));
    end
    f3=@(x)interp1(a,f1,x);
    X=fminsearch(@(x)abs(f3(x)-Level),0);
    res=f3(X*th/ThMax);
end