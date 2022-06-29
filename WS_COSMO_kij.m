function diff=WS_COSMO_kij(P,T,kij)
% required input parameters
% T temperature [K]; P pressure [Mpa]; component component parameters; kij
% output Ge difference
% adjust kij to give diff=0
c1={'74-98-6'};
c2={'75-46-7'};
load('database.mat');
for i=1:length(data)
    if isequal(data(i).cas,c1)
        break;
    end
end
for j=1:length(data)
    if isequal(data(j).cas,c2)
        break;
    end
end
component(1).Tc=data(i).tc;component(2).Tc=data(j).tc;
component(1).Pc=data(i).pc;component(2).Pc=data(j).pc;
component(1).Twu=data(i).c;component(2).Twu=data(j).c;
%component(1).omega=data(i).omega;component(2).omega=data(j).omega;
component(1).A=data(i).A;component(2).A=data(j).A;
component(1).V=data(i).V;component(2).V=data(j).V;
component(1).disp=data(i).disp;component(2).disp=data(j).disp;
component(1).sigma=data(i).sigma;component(2).sigma=data(j).sigma;
R=8.31433;
Tc=[component(1).Tc;component(2).Tc];
Pc=[component(1).Pc;component(2).Pc];
Twu=[component(1).Twu;component(2).Twu];
x=[0.5;0.5];
psigma_pure=[component(1).sigma;component(2).sigma]';
Area=[component(1).A;component(2).A];
Vol=[component(1).V;component(2).V];
ekB=[component(1).disp;component(2).disp];
alpha_i=(T./Tc).^(Twu(:,3).*(Twu(:,2)-1)).*exp(Twu(:,1).*(1-(T./Tc).^(Twu(:,3).*Twu(:,2))));
a_ci=0.45724*R*R*Tc.^2./Pc;
a_ii=alpha_i.*a_ci;
b_ii=0.0778*R*Tc./Pc;
b_aRT_ij=(repmat((b_ii-a_ii/R/T),[1,2])+repmat((b_ii-a_ii/R/T)',[2,1])).*[1,1-kij;1-kij,1]/2;
lngamma=GE_COSMO_SAC_dsp(T,x,Area,Vol,ekB,psigma_pure);
g_E_res=R*T*x'*lngamma;
Z_WS_COSMO_PURE=zeros(1,2);
%pure1
b1=b_ii(1);
a1=a_ii(1);
A1=a1*P/R/R/T/T;B1=b1*P/R/T;
C2=-(1-B1);C1=-(2*B1-A1+3*B1^2);C0=-(A1*B1-B1^2-B1^3);
root=roots([1,C2,C1,C0]);
if isreal(root(1)) && isreal(root(2)) && isreal(root(3))
Z_WS_COSMO_PURE(1)=min(root);
else
Z_WS_COSMO_PURE(1)=[isreal(root(1)),isreal(root(2)),isreal(root(3))]*root;
end
%pure2
b2=b_ii(2);
a2=a_ii(2);
A2=a2*P/R/R/T/T;B2=b2*P/R/T;
C2=-(1-B2);C1=-(2*B2-A2+3*B2^2);C0=-(A2*B2-B2^2-B2^3);
root=roots([1,C2,C1,C0]);
if isreal(root(1)) && isreal(root(2)) && isreal(root(3))
Z_WS_COSMO_PURE(2)=min(root);
else
Z_WS_COSMO_PURE(2)=[isreal(root(1)),isreal(root(2)),isreal(root(3))]*root;
end
%
b=x'*b_aRT_ij*x/(1-x'*(a_ii./b_ii/R/T)+g_E_res/R/T/0.62323);
a=b*(x'*(a_ii./b_ii)+(g_E_res)/-0.62323);
A=a*P/R/R/T/T;B=b*P/R/T;
C2=-(1-B);C1=-(2*B-A+3*B^2);C0=-(A*B-B^2-B^3);
root=roots([1,C2,C1,C0]);
if isreal(root(1)) && isreal(root(2)) && isreal(root(3))
Z_WS_COSMO=min(root);
else
Z_WS_COSMO=[isreal(root(1)),isreal(root(2)),isreal(root(3))]*root;
end
g_E_pr=Z_WS_COSMO-Z_WS_COSMO_PURE*x-log(Z_WS_COSMO./Z_WS_COSMO_PURE)*x-log((Z_WS_COSMO-B)./(Z_WS_COSMO_PURE-[B1,B2]))*x+A/2.82843/B*log((2*Z_WS_COSMO-0.82843*B)/(2*Z_WS_COSMO+4.82843*B))...
    -([A1,A2]/2.82843./[B1,B2].*log((2*Z_WS_COSMO_PURE-0.82843*[B1,B2])./(2*Z_WS_COSMO_PURE+4.82843*[B1,B2])))*x;
diff=g_E_res/R/T-g_E_pr;
end

