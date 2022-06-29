function lngamma=UNIFAC_PSRK(T,x,group1,group2)
% required input parameters
% T temperature; x concentration; group group lable [vk, k]
% output log(gamma)
%% load database
load UNIFAC_PSRK_BINARY.mat
load UNIFAC_PSRK_PURE.mat
%% get parameters
group_all=sort(unique([group1(:,2);group2(:,2)]));  
group1=[(group_all==group1(:,2)')*group1(:,1),group_all];
group2=[(group_all==group2(:,2)')*group2(:,1),group_all];
try
    R=(UNIFAC_PSRK_PURE(:,3)'*(UNIFAC_PSRK_PURE(:,1)==group_all'))';
    Q=(UNIFAC_PSRK_PURE(:,4)'*(UNIFAC_PSRK_PURE(:,1)==group_all'))';
    MAIN=(UNIFAC_PSRK_PURE(:,2)'*(UNIFAC_PSRK_PURE(:,1)==group_all'))';
catch
    lngamma=zeros(2,1);
    errordlg('group pure parameters missing','undefined');
    return
end

binary=zeros(length(group_all),length(group_all),3);
for i=1:length(group_all)
    for j=i+1:length(group_all)
        try
            label=UNIFAC_PSRK_BINARY(:,1)==MAIN(i) & UNIFAC_PSRK_BINARY(:,2)==MAIN(j);
            if sum(label)==0 && MAIN(i)~=MAIN(j)
               errordlg('group binary parameters missing','undefined');
               return
            else
                binary(i,j,[1,2,3])=UNIFAC_PSRK_BINARY(:,[3,4,5])'*label;
                binary(j,i,[1,2,3])=UNIFAC_PSRK_BINARY(:,[6,7,8])'*label;
            end
        catch
            lngamma=zeros(2,1);
            errordlg('group binary parameters missing','undefined');
            return
        end
    end
end
%% combinatorial contribution
ri=[group1(:,1)'*R;group2(:,1)'*R];
qi=[group1(:,1)'*Q;group2(:,1)'*Q];
phi_over_x=ri/(ri'*x);
theta_over_phi=qi./ri/(qi'*x)*(ri'*x);
li=10/2*(ri-qi)-(ri-1);
lngamma_C=log(phi_over_x)+10/2*qi.*log(theta_over_phi)+li-phi_over_x*(li'*x);
%% residual contribution
X_m=(group1(:,1)*x(1)+group2(:,1)*x(2))/(sum(group1(:,1))*x(1)+sum(group2(:,1))*x(2));
theta_m=Q.*X_m/(Q'*X_m);
theta_1=Q.*group1(:,1)/(Q'*group1(:,1));
theta_2=Q.*group2(:,1)/(Q'*group2(:,1));
phi=exp(-(binary(:,:,1)+binary(:,:,2)*T+binary(:,:,3)*T*T)/T);
%phi=[1 0.65 0.7366;1.193 1 0.8236; 1.05717 0.9829 1];
lngamma_group_m=Q.*(1-log(theta_m'*phi)'-(phi*(theta_m./(theta_m'*phi)')));
lngamma_group_1=Q.*(1-log(theta_1'*phi)'-(phi*(theta_1./(theta_1'*phi)')));
lngamma_group_2=Q.*(1-log(theta_2'*phi)'-(phi*(theta_2./(theta_2'*phi)')));
lngamma_R=zeros(2,1);
lngamma_R(1)=group1(:,1)'*(lngamma_group_m-lngamma_group_1);
lngamma_R(2)=group2(:,1)'*(lngamma_group_m-lngamma_group_2);
%% output
lngamma=lngamma_C+lngamma_R;
end

