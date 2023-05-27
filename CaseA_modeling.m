clear;
clc;
tstart=tic;

%%%%%%Parameter%%%%%%
%ep: small number;
%timeyear: coefficient convert year to second;
%g: gravitaional coefficient;
%rho: water density;
%R: submerged specific gravity;
%B: channel width;
%L: channel length;
%SI: initial slope;
%Se: equilibrium slope;
%Duration: duration of calculation (year);
%If: flood intermittency;
%Q: flow discharge (m3/s);
%qsfT: total volume sediment supply per unit width (m2/s);
%qsf: volume sediment supply per unit width for each size range (m2/s)
%Pbs_st: GSD of bed surface measured at St. Louis;
%Phi: phi value of grain size;
%D: characteristic grain size of each size range;
%ng: number of sediment group;
%Cf: resistance coefficient;
%na: coefficient relating active layer thickness;
%p: bed porosity;
%au: coefficient related to Exner discretization;
%Ls: thickness of storage layer;
%altr: coefficient related to interfacial exchange;
ep=1e-6;
If=0.73;  
timeyear=365.25*24*3600;
timeyear=timeyear*If;
g=9.81;
rho=1000;
R=1.65;
B=1300;    % width
L=1000000; % length
SI=0.0000489; 
Duration=18;  
Nt=18;
Q=26000;      % discharge

%upstream sediment inflow
Pbs_up=[52.3,43.3,4.4,0,0,0,0,0,0,0]/100;
Phin_up=[-3,-2,-1,0,1,2,3,4,5,6,7];
Phi_up=(Phin_up(1:10)+Phin_up(2:11))/2;
Phig_up=Pbs_up*Phi_up';
D=2.^Phi_up/1000;
ng=size(D,2);

Cf=0.00491;
p=0.35;
au=1;     
Ls=0.5;   
altr=0.5;  

%An et al.2021
a=0.037;
b=0.445;
c=0.876;
d=-0.348;

%dt: time step (second);
%dth: time step for hydraulic calculation (second);
%nh: dt/dth;
%dx: cell size;
%n: cell number;
%Courant: Courant number;
dt=0.001;
dth=0.00001;
nh=round(dt/dth);
n=101;
dx=L/(n-1);
Courant=zeros(n,1);
interval=round(Duration/dt/Nt);

Pbs_yz02=[27.9,70.3,1.8,0,0,0,0,0,0,0]/100;
Phin_yz02=[-3,-2,-1,0,1,2,3,4,5,6,7];
Phi_yz02=(Phin_yz02(1:10)+Phin_yz02(2:11))/2;
Phig_yz02=Pbs_yz02*Phi_yz02';
Dsg_02 = 2.^(Pbs_yz02*Phi_yz02')/1000;

Pbs_yz20=[4.40692,7.65467,17.43564,17.43564,8.44905,8.44905,8.44905,2.71999,22.23655,2.76345]/100;
Phin_yz20=[-3,-2,-1,0,1,2,3,4,5,6,7];
Phi_yz20=(Phin_yz20(1:10)+Phin_yz20(2:11))/2;
Phig_yz20=Pbs_yz20*Phi_yz20';


%drill sample, case A
Pbs_sub = [5.87590,10.20623,23.24752,23.24752,11.26540,11.26540,11.26540,3.62665,0,0]/100; 
Phin_yz = [-3,-2,-1,0,1,2,3,4,5,6,7];
Phi_yz = (Phin_yz(1:10)+Phin_yz(2:11))/2;
Dsg_sub = 2.^(Pbs_sub*Phi_yz')/1000;


%Pbs: GSD of sediment on bed surface;
%Dsg: geometrical mean grain size of sediment on bed surface
%Pbs=ones(n,1)*Pbs_st;
%Dsg=2.^(Pbs*Phi')/1000;
%qw: flow discharge per unit width;
%h: water depth;
%u: flow velocity;
%two ghost cells are speficied at both upstream boundary and downstream boundary;
qw=ones(n+2,1)*Q/B;
h=(Cf.*qw.^2./g/SI).^(1/3);
u=qw(2:n+1)./h(2:n+1);
%taub: bed shear stress;
%ustar: shear velocity;
%shi: shields number for sediment of each size range;
%NA: coefficient in Naito's relation;
%NB: exponent in Naito's relation;
%Nstar: dimensionless sediment transport rate;
%qse: equilibrium sediment transport rate per unit width for each size range;
%qs: sediment transport rate per unit width;
%qsT: total volume sediment transport rate per unit width;m2/s
%Gta: ambient annual sediment load;
taub=rho.*Cf.*u.^2;
ustar=u.*(Cf)^0.5;
shi=taub*ones(1,ng)./rho./R./g./(ones(n,1)*D);
NA=a.*((ones(n,1)*D)./(Dsg_02*ones(1,ng))).^b;
NB=c.*((ones(n,1)*D)./(Dsg_02*ones(1,ng))).^d;
Nstar=NA.*(shi).^NB;
%Pbs_tmp = Pbs_yz02(:,1:3);
qse=Nstar.*Pbs_yz02.*((ustar.^3)*ones(1,ng))./R./g./Cf;
%qse=qse*0.98;
qs=qse;          
qsT=sum(qs,2);  
Pstr=qs./(qsT*ones(1,ng));  
Dg_load=2.^(Pstr*Phi_up')/1000;
Gta=qsT/1000*timeyear*rho*(R+1)*B/1e6;
%Gta(1e6 ton)=qsT(m2/s)*t(Ò»Äês)*2650(kg/m3)*B(m)/1e9
%qsfT_bench: benchmark total volume sediment supply per unit width (m2/s);
%qsf_bench: benchmark volume sediment supply per unit width for each size range (m2/s)
%qsfT_bench=qsT(1);
%qsf_bench=qs(1,:);


%zb: bed elevation;
%zbini: initial bed elevation;
%s: bed slope;
%FR: Froude number;
x=linspace(0,L,n)';
zb=SI*(L-x);
zbini=zb;
s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
h=(Cf.*qw.^2./g/SI).^(1/3);
u=qw(2:n+1)./h(2:n+1);
FR=u./(g*h(2:n+1)).^0.5;
%La: thickness of active layer;
%Psub: GSD of substrate sediment;
%Store: information of substrate stratigraphy;
%indup: number of sublayer at every node;
%Pup:proportion of sediment on uppermost sublayer;
%Lup: thickness of uppermost sublayer
La=ones(n,1)*1.5;
Psub=Pbs_sub;
Store=cell(n,1);
indup=ceil((zb-La+10*Ls)./Ls);
for is=1:n 
    Store{is,1}=zeros(100,ng);
    Store{is,1}(1:indup(is),:)=ones(indup(is),1)*Psub;
end
Pup=ones(n,1)*Psub;
Lup=zb+10*Ls-La-Ls*(indup-1);
PLa=Pbs_yz02;
Store_ini=Store;
t=0;
i=0;


%left&right conditions for flux calculation
hl=zeros(n+1,1);
hr=zeros(n+1,1);
qwl=zeros(n+1,1);
qwr=zeros(n+1,1);
ul=zeros(n+1,1);
ur=zeros(n+1,1);
%numerical flux
Fl=zeros(n+1,2);
Fr=zeros(n+1,2);
Fhll=zeros(n+1,2);
F=zeros(n+1,2);
%source terms
ss=zeros(n,1);
sf=zeros(n,1);
%fI: interfacial exchange fractions;
%dzb: change of bed evolution;
%dLa: change of activer layer thickness;
%dPbs: change of surface fraction;
%delta: change of substrate elevation;
fI=zeros(n,ng);
dzb=zeros(n,1);
dLa=zeros(n,1);
dPLa=zeros(n,ng);
delta=zeros(n,1);
%qwt: flow discharge per unit width at different time;
%zbt: bed elevation at the end of each year;
%st: bed slope at the end of each year;
%Dsgt: Dsg at different time;
%Pbst: Pbs at different time;
%qsTt: total sediment transport rate per unit width at different time;
%Gtt: annual sediment load at different time;
%Dg_loadt: Dg_load at different time;
qwt=zeros(n,Nt);
ht=zeros(n,Nt);
zbt=zeros(n,Nt);
st=zeros(n,Nt);
Dsgt=zeros(n,Nt);
%Pbst=cell(Nt,1);
qsTt=zeros(n,Nt);
qst1=zeros(n,Nt);
qst2=zeros(n,Nt);
qst3=zeros(n,Nt);
qst4=zeros(n,Nt);
qst5=zeros(n,Nt);
qst6=zeros(n,Nt);
Gtt=zeros(n,Nt);
Nstar_total_t=zeros(n,Nt);
Dg_loadt=zeros(n,Nt);
it=1;
F(1,1)=qw(1);

%%%%%%Time Processing%%%%%%
while t<Duration
    %%%%%%Flow Hydraulics£ºShallow Water Equation%%%%%%
    for ih=1:nh
        ih;
        %%%BC: upstream given discharge and mass balance, downstream bed resistance%%%
        qw(1)=Q/B;
        hsub=h(1)+(qw(1)-F(1,1))*dth*timeyear/dx;
        hsuper=(Cf.*qw(1).^2./g/s(1)).^(1/3);
        h(1)=(FR(1)>=1)*hsuper+(FR(1)<1)*hsub;
        qw(n+2)=qw(n+1);
        h(n+2)=(Cf.*qw(n+2).^2./g/s(n)).^(1/3);
        %%%HLL%%%
        %left&right conditions
        hl=h(1:n+1);
        hr=h(2:n+2);
        qwl=qw(1:n+1);
        qwr=qw(2:n+2);
        ul=qwl./hl;
        ur=qwr./hr;
        %wave speed
        cal=sqrt(g*hl);
        car=sqrt(g*hr);
        us=0.5*(ul+ur)+cal-car;
        cas=0.5*(cal+car)+0.25*(ul-ur);
        sl=min(ul-cal,us-cas)*ones(1,2);
        sr=max(ur+car,us+cas)*ones(1,2);
        %numerical flux
        Fl=[qwl,qwl.^2./hl+0.5*g*hl.^2];
        Fr=[qwr,qwr.^2./hr+0.5*g*hr.^2];
        Fhll=(sr.*Fl-sl.*Fr+sl.*sr.*[hr-hl,qwr-qwl])./(sr-sl);
        F=(sl<0).*(sr>0).*Fhll+(sl>=0).*Fl+(sr<=0).*Fr;
        %%%source term%%
        ss=g.*h(2:n+1).*s;
        sf=Cf*abs(qw(2:n+1))./h(2:n+1).^2;
        %%%time advance%%%
        h(2:n+1)=h(2:n+1)+(F(1:n,1)-F(2:n+1,1))/dx*dth*timeyear;
        qw(2:n+1)=(qw(2:n+1)+(F(1:n,2)-F(2:n+1,2))/dx*dth*timeyear+ss*dth*timeyear)./(1+sf*dth*timeyear);
        Courant=max(Courant,(qw(2:n+1)./h(2:n+1)+(g.*h(2:n+1)).^0.5)*dth*timeyear/dx);
    end
    u=qw(2:n+1)./h(2:n+1);
    FR=u./(g*h(2:n+1)).^0.5;
    

    %%%sediment transport
    taub=rho.*Cf.*u.^2;
    ustar=u.*(Cf)^0.5;
    shi=taub*ones(1,ng)./rho./R./g./(ones(n,1)*D);
    NA=a.*((ones(n,1)*D)./(Dsg_02*ones(1,ng))).^b;
    NB=c.*((ones(n,1)*D)./(Dsg_02*ones(1,ng))).^d;
    Nstar=NA.*(shi).^NB;
    Nstar_total = sum(Nstar,2);
    %Pbs_tmp = Pbs(1:3);
    qse=Nstar.*PLa.*((ustar.^3)*ones(1,ng))./R./g./Cf;
    %qse=qse*0.98;
    qs=qse;
    qsT=sum(qs,2);
    Pstr=qs./(qsT*ones(1,ng));
    Dg_load=2.^(Pstr*Phi_up')/1000;
    Gt=qsT/1000*timeyear*rho*(R+1)*B/1e6;
      
 

    qsfT=4.82*10^(-5);  
    qsf=qsfT*Pbs_up;  
    qsTback=[qsfT;qsT(1:n-1)];
    qsTit=qsT;
    qsTfrnt=[qsT(2:n);2*qsT(n)-qsT(n-1)];
    qsTdif=au*(qsTit-qsTback)+(1-au)*(qsTfrnt-qsTit);
    dzb=-qsTdif/dx*dt*timeyear/(1-p);
    dzb(n)=0;
    zb=zb+dzb;
    %%%evolution of surface fraction
    delta=dzb-dLa;
    %zeros_tmp = zeros(101, 6);
    %Pstr_tmp = [Pstr zeros_tmp];
    fI=((delta<=0)*ones(1,ng)).*Pup+((delta>0)*ones(1,ng)).*(altr*PLa+(1-altr)*Pstr);
    qsback=[qsf;qs(1:n-1,:)];
    qsit=qs;
    qsfrnt=[qs(2:n,:);2*qs(n,:)-qs(n-1,:)];
    qsdif=au*(qsit-qsback)+(1-au)*(qsfrnt-qsit);
    %qsdif_tmp = [qsdif zeros_tmp];
    dPLa=((-qsdif+qsTdif*ones(1,ng).*fI)/dx*dt*timeyear/(1-p)-(PLa-fI).*(dLa*ones(1,ng)))./(La*ones(1,ng));
    %Pbs;
    PLa=PLa+dPLa;
    PLa=(PLa>0).*PLa;
    PLa=PLa./(sum(PLa,2)*ones(1,ng));
    PLa;
    
    %%%Stratigraphy storage
    indn=find((delta<=Ls-Lup)&(delta>=-Lup));
    indinc=find(delta>Ls-Lup);
    inddec=find(delta<-Lup);
    Pup(indn,:)=(Pup(indn,:).*(Lup(indn)*ones(1,ng))+fI(indn,:).*(delta(indn)*ones(1,ng)))./((Lup(indn)+delta(indn))*ones(1,ng));
    %Pup=(Pup>0).*Pup;
    %Pup=Pup./(sum(Pup,2)*ones(1,ng));
    Lup(indn)=Lup(indn)+delta(indn);
    if size(indinc,1)>0
        for is=1:size(indinc,1)
            ii=indinc(is);
            Pup(ii,:)=(Pup(ii,:).*Lup(ii)+fI(ii,:).*(Ls-Lup(ii)))./Ls;
            %Pup(ii,:)=(Pup(ii,:)>0).*Pup(ii,:);
            %Pup(ii,:)=Pup(ii,:)./sum(Pup(ii,:),2);
            Store{ii,1}(indup(ii),:)=Pup(ii,:);
            inc=ceil((delta(ii)-(Ls-Lup(ii)))/Ls);
            Store{ii,1}(indup(ii)+1:indup(ii)+inc,:)=ones(inc,1)*fI(ii,:);
            indup(ii)=indup(ii)+inc;
            Pup(ii,:)=fI(ii,:);
            Lup(ii)=delta(ii)-(Ls-Lup(ii))-(inc-1)*Ls;
        end
    end
    if size(inddec,1)>0
        for is=1:size(inddec,1)
            id=inddec(is);
            dec=ceil((-Lup(id)-delta(id))/Ls);
            Store{id,1}(indup(id)-dec+1:indup(id),:)=zeros(dec,ng);
            indup(id)=indup(id)-dec;
            Pup(id,:)=Store{id,1}(indup(id),:);
            Lup(id)=dec*Ls+Lup(id)+delta(id);
        end
    end
    %%%parameter update
    t=t+dt;
    i=i+1;
    s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
    Dsg_02=2.^(PLa*Phi_yz')/1000;
    
    %%%record information evry several years
    if (mod(i,interval)==0)
        PLa;
        qwt(:,it)=qw(2:n+1);
        ht(:,it)=h(2:n+1);
        zbt(:,it)=zb;
        st(:,it)=s;
        Dsgt(:,it)=Dsg_02;
        %Pbst{it,1}=Pbs;
        qsTt(:,it)=qsT;
        qst1(:,it)=qs(:,1);
        qst2(:,it)=qs(:,2);
        qst3(:,it)=qs(:,3);
        qst4(:,it)=qs(:,4);
        qst5(:,it)=qs(:,5);
        qst6(:,it)=qs(:,6);
        Gtt(:,it)=Gt;
        Nstar_total_t(:,it)=Nstar_total;
        Dg_loadt(:,it)=Dg_load;
        it
        it=it+1;
    end
    
    %figure(1);plot(x,qw(2:n+1),'*r')
    %figure(1);plot(x,h(2:n+1),'*r')
    %figure(1);plot(x,zb,'*r')
    %figure(1);plot(x,sc,'*r')
    %figure(1);plot(x,Dsg,'*r')
    %figure(1);plot(x,Dg_load,'*r')
    %figure(1);plot(x,Gt,'*r)
    %drawnow
    
end
PLa

for is=1:n
    Store{is,1}(indup(is),:)=Pup(is,:);
end
tend=toc(tstart);
tc=tend/60