function ln_phi=WS_UNIFAC_PR(P,T,component,phase_flag,kij)
% required input parameters
% T temperature [K]; P pressure [MPa]; component component parameters; phase_flag
% vapor(1)/liquid(0)
% output log(phi)

R=8.31433;
Tc=[component(1).Tc;component(2).Tc];
Pc=[component(1).Pc;component(2).Pc];
Twu=[component(1).Twu;component(2).Twu];
x=(1-phase_flag)*[component(1).x;component(2).x]+phase_flag*[component(1).y;component(2).y];
group1=component(1).group;
group2=component(2).group;
Tr=T./Tc;
alpha_i=(T./Tc).^(Twu(:,3).*(Twu(:,2)-1)).*exp(Twu(:,1).*(1-(T./Tc).^(Twu(:,3).*Twu(:,2))));
a_ci=0.45724*R*R*Tc.^2./Pc;
a_ii=alpha_i.*a_ci;
b_ii=0.0778*R*Tc./Pc;
b_aRT_ij=(repmat((b_ii-a_ii/R/T),[1,2])+repmat((b_ii-a_ii/R/T)',[2,1])).*[1,1-kij;1-kij,1]/2;
lngamma=UNIFAC_DO(T,x,group1,group2);
g_E_res=R*T*x'*lngamma;
b=x'*b_aRT_ij*x/(1-x'*(a_ii./b_ii/R/T)+g_E_res/R/T/0.62323);
a=b*(x'*(a_ii./b_ii)+(g_E_res)/-0.62323);
A=a*P/R/R/T/T;B=b*P/R/T;
C2=-(1-B);C1=-(2*B-A+3*B^2);C0=-(A*B-B^2-B^3);
root=roots([1,C2,C1,C0]);
if isreal(root(1)) && isreal(root(2)) && isreal(root(3))
Z_WS_UNI=phase_flag*max(root)+(1-phase_flag)*min(root);
else
Z_WS_UNI=[isreal(root(1)),isreal(root(2)),isreal(root(3))]*root;
end
Z_real=Z_WS_UNI;
QQ=x'*b_aRT_ij*x;
DD=x'*(a_ii./b_ii/R/T)-g_E_res/R/T/0.62323;
dQQ=2*b_aRT_ij*x;
dDD=a_ii./b_ii/R/T+lngamma/-0.62323;
db=1/(1-DD).*dQQ-QQ/(1-DD)/(1-DD)*(1-dDD);
da=R*T*DD*db+R*T*b*dDD;
ln_phi_0=(db/b)*(Z_real-1);
ln_phi_1=log(Z_real-B);
ln_phi_2=da/a-db/b;
ln_phi_3=A/B/2.82843*log((2*Z_real-0.82843*B)/(2*Z_real+4.82843*B));
ln_phi=ln_phi_0-ln_phi_1+ln_phi_2*ln_phi_3;
end

