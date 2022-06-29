function ln_phi=PSRK(P,T,component,phase_flag)
% required input parameters
% T temperature [K]; P pressure [MPa]; component component parameters; phase_flag
% vapor(1)/liquid(0)
% output log(phi)

R=8.31433;
Tc=[component(1).Tc;component(2).Tc];
Pc=[component(1).Pc;component(2).Pc];
omega=[component(1).omega;component(2).omega];
x=(1-phase_flag)*[component(1).x;component(2).x]+phase_flag*[component(1).y;component(2).y];
group1=component(1).group;
group2=component(2).group;
Tr=T./Tc;
alpha_i=(1+(0.48+1.574*omega-0.176*omega.^2).*(1-Tr.^0.5)).^2;
a_ci=0.42748*R*R*Tc.^2./Pc;
a_ii=alpha_i.*a_ci;
b_ii=0.08664*R*Tc./Pc;
b=x'*b_ii;
lngamma=UNIFAC_PSRK(T,x,group1,group2);
g_E_res=R*T*x'*lngamma;
a=b*(x'*(a_ii./b_ii)+(g_E_res+R*T*x'*log(b./b_ii))/-0.64663);
A=a*P/R/R/T/T;B=b*P/R/T;
C2=-1;C1=A-B-B^2;C0=-A*B;
root=roots([1,C2,C1,C0]);
if isreal(root(1)) && isreal(root(2)) && isreal(root(3))
Z_SRK=phase_flag*max(root)+(1-phase_flag)*min(root);
else
Z_SRK=[isreal(root(1)),isreal(root(2)),isreal(root(3))]*root;
end
Z_real=Z_SRK;
ln_phi_0=(b_ii/b)*(Z_real-1);
ln_phi_1=log(Z_real-B);
ln_phi_2=a_ii./b_ii/R/T+lngamma/-0.64663+(log(b./b_ii)+b_ii/b-1)/-0.64663;
ln_phi_3=log((Z_real)/(Z_real+B));
ln_phi=ln_phi_0-ln_phi_1+ln_phi_2*ln_phi_3;
end

