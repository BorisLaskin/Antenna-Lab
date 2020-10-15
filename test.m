J0=@(x)1/pi*integral(@(t)cos(x.*sin(t)),0,pi);
a=linspace(0,100,1000);
f1=zeros(size(a));
for i=1:numel(a)
    f1(i)=J0(a(i));
end
f2=@(x)interp1(a,f1.^2,x);
f3=@(x)interp1(a,f1,x);
Full=integral(f2,0,100);
X=@(Level)fminsearch(@(x)abs(f3(x)-Level),0);
Part=@(Level)(integral(@(x)f3(x),0,X(Level))).^2;
Fres=@(L)Part(L)./Full;

E1=@(D,L)(abs(L.^2 - 1).*abs(D)^2)./(4*abs(acos(L)).^2);
E2=@(D,L)D/4 + (D*L.*(1 - L.^2).^(1/2))./(4*acos(L));
E3=@(D,L)2*E1(D,L)./E2(D,L)/D;

Res=@(D,L)E3(D,L).*Fres(L);
L1=linspace(0,1,100);
Result=zeros(size(L1));
for i=1:numel(L1)
    Result(i)=Res(10,L1(i));
end
plot(L1,Result); % оптимальное значение level вне зависимости от диаметра равняется 0,17