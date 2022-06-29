function lngamma=GE_COSMO_SAC_dsp(T,x,A,V,ekB,psigma_pure)
% required input parameters
% T temperature; x concentration; A surface area; V volume; 
% ekB dispersion parameters psigma_pure sigma profiles 
% output log(gamma)
a_eff=7.25;
psigma_pure=psigma_pure./sum(psigma_pure);
lngamma_comb_i=lngamma_comb(x,A,V);
lngamma_disp_i=lngamma_disp(x,ekB);
psigma_s=psigma_pure*(x.*A)/(x'*A);
lngamma_s=lngamma_resid(T,psigma_s);
lngamma_1=lngamma_resid(T,psigma_pure(:,1));
lngamma_2=lngamma_resid(T,psigma_pure(:,2));
lngamma_resid_1=A(1)/a_eff*(psigma_pure(:,1)'*(lngamma_s-lngamma_1));
lngamma_resid_2=A(2)/a_eff*(psigma_pure(:,2)'*(lngamma_s-lngamma_2));
lngamma=lngamma_comb_i+lngamma_disp_i+[lngamma_resid_1;lngamma_resid_2];
end
%% the combinatorial part of ln(gamma_i)
function y=lngamma_comb(x,A,V) % input x concentration A surface area V volume
q0=79.53; %normalized surface area parameter
r0=66.69; %normalized volume parameter
z_coordination=10;     %coordination number
q=A/q0; r=V/r0;
l=z_coordination/2*(r-q)-(r-1);
phi_i_over_x_i=r/(x'*r);
theta_i_over_phi_i=q./r/(x'*q)*(x'*r);
y=log(phi_i_over_x_i)+z_coordination/2*q.*log(theta_i_over_phi_i)+l-phi_i_over_x_i*(x'*l);
end
%% the dispersion part of ln(gamma_i)
function y=lngamma_disp(x,ekB) % input x concentration ekB dispersion parameters
y=zeros(2,1);% initilize
w=0.27027; %default value
A=w*(0.5*sum(ekB)-(ekB(1)*ekB(2))^0.5);
y(1)=A*x(2)*x(2);y(2)=A*x(1)*x(1);
end
%% the residual part of ln(gamma_i)
function y=lngamma_resid(T,psigma) %input T temperature psigma sigma profiles
gamma_old=ones(length(psigma),1);%initialization
%constant
A_ES=6525.69; B_ES=1.4859*10^8; c_OT_OT=932.31;R=8.314/4184; max_iter=2000;counter=0;
C_ES=A_ES+B_ES/T/T;

if length(psigma)==51
    sigma=(-0.025:0.001:0.025)';
    Delta_W=-(C_ES*(repmat(sigma,1,51)+repmat(sigma',51,1)).^2)/R/T;
elseif length(psigma)==102
    sigma=[(-0.025:0.001:0.025),(-0.025:0.001:0.025)]';
    chb=zeros(102,102);chb(52:end,52:end)=[zeros(25,26),ones(25,25);zeros(1,51);ones(25,25),zeros(25,26)];
    Delta_W=-(C_ES*(repmat(sigma(:,1),1,102)+repmat(sigma(:,1)',102,1)).^2-c_OT_OT*(chb.*(repmat(sigma(:,1),1,102)-repmat(sigma(:,1)',102,1)).^2))/R/T;
end
gamma_new=((exp(Delta_W)*(psigma.*gamma_old)).^-1+gamma_old)/2;
while max(abs(gamma_new-gamma_old))>1e-7 && counter<max_iter
    counter=counter+1;
    gamma_old=gamma_new;
    gamma_new=((exp(Delta_W)*(psigma.*gamma_old)).^-1+gamma_old)/2;
end
if counter >=max_iter 
    errordlg('max step reached','Iteration error')
    y=zeros(length(psigma),1);
else 
    y=log(gamma_new);
end
end
