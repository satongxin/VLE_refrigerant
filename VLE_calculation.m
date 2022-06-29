function [y_cal,P_cal]=VLE_calculation
% this is the main function to perform the VLE calculation
% specify the CAS number of your binary mixture
c1="76-16-4";
c2="420-46-2";
% provide experimental data for initialization and error evaluation
expe=[0	0.065	0.1222	0.1839	0.2606	0.3672	0.454	0.5543	0.7406	0.8043	0.8754	0.9045	1
273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3	273.3
0	0.2392	0.3647	0.4568	0.5364	0.6171	0.6678	0.7209	0.8153	0.8507	0.8966	0.9168	1
632500	811900	947400	1075200	1210600	1357200	1463400	1572900	1738300	1782400	1824300	1838000	1850400
];
x_exp=expe(1,:); %x in mole fraction
T_exp=expe(2,:); %T in K
y_exp=expe(3,:); %y in moel fraction
P_exp=expe(4,:)/1000000; % P in MPa
%% look for pure component parameters
load('database.mat'); % load predefined parameter database
for i=1:length(data)      %match component in the database according to the CAS number
    if isequal(string(data(i).cas),c1)
        break;
    end
end
for j=1:length(data)
    if isequal(string(data(j).cas),c2)
        break;
    end
end
% read parameters
component(1).Tc=data(i).tc;component(2).Tc=data(j).tc;    %Tc
component(1).Pc=data(i).pc;component(2).Pc=data(j).pc;    %Pc
component(1).Twu=data(i).c;component(2).Twu=data(j).c;    %Twu parameters
component(1).omega=data(i).omega;component(2).omega=data(j).omega;   %omega
component(1).group=data(i).group1;component(2).group=data(j).group1;  %group-VTPR group1-PSRK group2-UNIFAC
component(1).A=data(i).A;component(2).A=data(j).A;  %area for COSMO calculation
component(1).V=data(i).V;component(2).V=data(j).V;  %volume for COSMO calculation
component(1).disp=data(i).disp; %dispersion parameter for COSMO calculation
component(2).disp=data(j).disp; %dispersion parameter for COSMO calculation
component(1).sigma=data(i).sigma;component(2).sigma=data(j).sigma; %sigma profile
%% iterative bubble point calculation
len=length(x_exp); 
y_cal=zeros(size(y_exp));
P_cal=zeros(size(P_exp));
for i=1:len
    P=P_exp(i);
    T=T_exp(i);
    component(1).x=x_exp(i);component(2).x=1-x_exp(i);
    component(1).y=y_exp(i);component(2).y=1-y_exp(i);
    eps=1/10^6;  
    dlp=500/10^6;
    k=1;
    while abs(dlp)>eps  %stop if less than tolerance
        k=k+1;
        if k>400 %stop if reach the max step
            y(1)=0; 
            break;
        end
        s1=sumy(P,T,component);
        P=P-eps;
        s2=sumy(P,T,component);
        dsdp=(sum(s2)-sum(s1))/eps;
        dlp=(sum(s2)-1)/dsdp;
        if dlp>0.1*P  
            dlp=0.1*P;
        end
        if dlp<-0.1*P
            dlp=-0.1*P;
        end
        P=P+dlp+eps/2;
        if P>2*P_exp(i)
            P=2*P_exp(i);
        end
        if P<0.2*P_exp(i)
            P=0.2*P_exp(i);
        end
        y=s2/sum(s2);
        component(1).y=y(1);component(2).y=y(2);
    end    
    y_cal(i)=y(1);P_cal(i)=P;
end
%compare calculated and experimental values in a plot
plot(x_exp,P_exp,'bo',y_exp,P_exp,'bo');hold on;plot(x_exp,P_cal,'k-',y_cal,P_cal,'k-');
end
function y_new=sumy(P,T,component)
% switch on one of the thermodynamic methods below
%the kij parameter in the WS mixing rule is calculated using 
%"WS_COSMO_kij" or "WS_UNIFAC_kij" 

%ln_phi_x=VTPR(P,T,component,0);
%ln_phi_y=VTPR(P,T,component,1);
%ln_phi_x=PSRK(P,T,component,0);
%ln_phi_y=PSRK(P,T,component,1);
%ln_phi_x=HVO(P,T,component,0);
%ln_phi_y=HVO(P,T,component,1);
ln_phi_x=MHV1(P,T,component,0);
ln_phi_y=MHV1(P,T,component,1);
%ln_phi_x=WS_UNIFAC_PR(P,T,component,0,0.007);
%ln_phi_y=WS_UNIFAC_PR(P,T,component,1,0.007);
%ln_phi_x=WS_COSMO_PR(P,T,component,0,0.2);
%ln_phi_y=WS_COSMO_PR(P,T,component,1,0.2);
y_new=[component(1).x;component(2).x].*exp(ln_phi_x-ln_phi_y);
end