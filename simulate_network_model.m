function [] = simulate_network_model(IT,pd,corstim,pick_dbs_freq)

%IT - iteration number (trial no)
%pd - 0(normal/healthy condition), 1(Parkinson's disease(PD) condition)
%corstim (cortical stimulation) - 0(off), 1(on)
%pick_dbs_freq - choose appropriate DBS frequency

rng shuffle

n = 10;             % number of neurons in each nucleus
tmax = 2000;       % ms
dt = 0.01;          % ms
t=0:dt:tmax;


% DBS Parameters

PW = 0.3;           % ms [DBS pulse width]
amplitude = 300;    % nA/cm2 [DBS current amplitude]
freqs=0:5:200;      % DBS frequency in Hz
pattern = freqs(pick_dbs_freq);

% Create DBS Current

if pick_dbs_freq==1
    
  Idbs=zeros(1,length(t)); 
  
else

  Idbs=creatdbs(pattern,tmax,dt,PW,amplitude);
  
end

% Create Cortical Stimulus Pulse

if corstim==1
    
  Iappco=zeros(1,length(t));
  Iappco((1000/dt):((1000+0.3)/dt))=350;
  
else

  Iappco=zeros(1,length(t));
  
end

% Run CTX-BG-TH Network Model
[TH_APs STN_APs GPe_APs GPi_APs Striat_APs_indr Striat_APs_dr Cor_APs] = CTX_BG_TH_network(pd,corstim,tmax,dt,n,Idbs,Iappco);

% Calculate GPi pathological low-frequency oscillatory power
dt1=0.01*10^-3;
params.Fs = 1/dt1; %Hz
params.fpass = [1 100];
params.tapers = [3 5];
params.trialave = 1;

[gpi_alpha_beta_area gpi_S gpi_f]=make_Spectrum(GPi_APs,params);

% gpi_alpha_beta_area - GPi spectral power integrated in 7-35Hz band
% gpi_S - GPi spectral power
% gpi_f - GPi spectral frequencies
 

 name = [num2str(IT) 'pd' num2str(pattern) '.mat'];
 eval(['save ' name])



quit

end

function Idbs=creatdbs(pattern,tmax,dt,PW,amplitude)

t=0:dt:tmax; Idbs=zeros(1,length(t)); 
iD=amplitude;
pulse=iD*ones(1,PW/dt);



i=1;
while i<length(t)
      
    Idbs(i:i+PW/dt-1)=pulse;
    instfreq=pattern;
    isi=1000/instfreq;
    i=i+round(isi*1/dt);
 end

end



function [TH_APs STN_APs GPe_APs GPi_APs Striat_APs_indr Striat_APs_dr Cor_APs] = CTX_BG_TH_network(pd,corstim,tmax,dt,n,Idbs,Iappco)
   
    %time step
    t=0:dt:tmax;

    %initial conditions (IC) to the different cell types
    v1=-62+randn(n,1)*5;    % normal distribution with mean = -62/-63.8 and SD = 5
    v2=-62+randn(n,1)*5;
    v3=-62+randn(n,1)*5;
    v4=-62+randn(n,1)*5;
    v5=-63.8+randn(n,1)*5;
    v6=-63.8+randn(n,1)*5;

% IC-Excitatory Regular Spiking Neurons
ae=0.02;
be=0.2;
ce=-65;
de=8;

% IC-Inhibitory Fast Spiking InterNeurons
ai=0.1;
bi=0.2;
ci=-65;
di=2;



Cm=1; %Membrane capacitance


%Ionic conductance and Equilibrium potential values
gl=[0.05 0.35 0.1 0.1]; El=[-70 -60 -65 -67];
gna=[3 49 120 100]; Ena=[50 60 55 50]; 
gk=[5 57 30 80]; Ek=[-75 -90 -80 -100];
gt=[5 5 0.5]; Et=0;
gca=[0 2 0.15]; Eca=[0 140 120];
Em=-100; 
gahp=[0 20 10]; k1=[0 15 10]; kca=[0 22.5 15];
ga=5;
gL=15;
gcak=1;

Kca=2*10^-3;
Z=2;
F=96485;
CAsn2=0.005*ones(n,1);
Cao=2000;
R=8314;
T=298;

alp=1/(Z*F);
con=(R*T)/(Z*F);




%%Setting initial matrices
vth=zeros(n,length(t)); %thalamic membrane voltage
vsn=zeros(n,length(t)); %STN membrane voltage
vge=zeros(n,length(t)); %GPe membrane voltage
vgi=zeros(n,length(t)); %GPi membrane voltage
vstr_indr=zeros(n,length(t)); %Indirect Striatum membrane voltage
vstr_dr=zeros(n,length(t)); %Direct Striatum membrane voltage
ve=zeros(n,length(t)); %Excitatory Cortex membrane voltage
vi=zeros(n,length(t)); %Inhibitory Cortex membrane voltage
ue=zeros(n,length(t));
ui=zeros(n,length(t));



%%initial conditions
vth(:,1)=v1;
vsn(:,1)=v2;
vge(:,1)=v3;
vgi(:,1)=v4;
vstr_indr(:,1)=v5;
vstr_dr(:,1)=v6;
ve(:,1)=ce;
ue(:,1)=be*ve(1);
vi(:,1)=ci;
ui(:,1)=bi*vi(1);

% State variables 
N3=gpe_ninf(vge(:,1));N4=gpe_ninf(vgi(:,1));
H1=th_hinf(vth(:,1)); 
H3=gpe_hinf(vge(:,1));H4=gpe_hinf(vgi(:,1));
R1=th_rinf(vth(:,1)); 
R3=gpe_rinf(vge(:,1));R4=gpe_rinf(vgi(:,1));
CA2=0.1; 
CA3=CA2;CA4=CA2; 

N2=stn_ninf(vsn(:,1));
H2=stn_hinf(vsn(:,1));
M2=stn_minf(vsn(:,1));
A2=stn_ainf(vsn(:,1));
B2=stn_binf(vsn(:,1));
C2=stn_cinf(vsn(:,1));
D2=stn_d2inf(vsn(:,1));
D1=stn_d1inf(vsn(:,1));
P2=stn_pinf(vsn(:,1));
Q2=stn_qinf(vsn(:,1));
R2=stn_rinf(vsn(:,1));

m5=alpham(vstr_indr(:,1))./(alpham(vstr_indr(:,1))+betam(vstr_indr(:,1)));
h5=alphah(vstr_indr(:,1))./(alphah(vstr_indr(:,1))+betah(vstr_indr(:,1)));
n5=alphan(vstr_indr(:,1))./(alphan(vstr_indr(:,1))+betan(vstr_indr(:,1)));
p5=alphap(vstr_indr(:,1))./(alphap(vstr_indr(:,1))+betap(vstr_indr(:,1)));
m6=alpham(vstr_dr(:,1))./(alpham(vstr_dr(:,1))+betam(vstr_dr(:,1)));
h6=alphah(vstr_dr(:,1))./(alphah(vstr_dr(:,1))+betah(vstr_dr(:,1)));
n6=alphan(vstr_dr(:,1))./(alphan(vstr_dr(:,1))+betan(vstr_dr(:,1)));
p6=alphap(vstr_dr(:,1))./(alphap(vstr_dr(:,1))+betap(vstr_dr(:,1))); 

%%Synapse parameters
Esyn = [-85 0 -85 0 -85 -85 -80];
tau=5; tau_i=13; gpeak=0.43; gpeak1=0.3; 

S2a=zeros(n,1); 
S21a=zeros(n,1); 
S2b=zeros(n,1); 
S21b=zeros(n,1); 
S2an=zeros(n,1); 
S21an=zeros(n,1); 
S3a=zeros(n,1); 
S31a=zeros(n,1); 
S3b=zeros(n,1);
S31b=zeros(n,1);
S32b=zeros(n,1);
S3c=zeros(n,1);
S31c=zeros(n,1); 
S32c=zeros(n,1); 
S4=zeros(n,1); 
S5=zeros(n,1);
S51=zeros(n,1);
S52=zeros(n,1);
S53=zeros(n,1);
S54=zeros(n,1);
S55=zeros(n,1);
S56=zeros(n,1);
S57=zeros(n,1);
S58=zeros(n,1);
S59=zeros(n,1);
S9=zeros(n,1);
S6a=zeros(n,1);
S6b=zeros(n,1);
S6bn=zeros(n,1);
S61bn=zeros(n,1);
S61b=zeros(n,1);
S91=zeros(n,1);
S92=zeros(n,1);
S93=zeros(n,1);
S94=zeros(n,1);
S95=zeros(n,1);
S96=zeros(n,1);
S97=zeros(n,1);
S98=zeros(n,1);
S99=zeros(n,1);
S7=zeros(n,1);
S8=zeros(n,1);
S1a=zeros(n,1); 
S1b=zeros(n,1); 
S1c=zeros(n,1); 
Z1a=zeros(n,1);
Z1b=zeros(n,1);

t_a = 1000; % Max duration of syn conductance
t_vec = 0:dt:t_a;
const = gpeak/(tau*exp(-1));
const1 = gpeak1/(tau*exp(-1)); 
const2 = gpeak1/(tau*exp(-1)); 

%TH-CTX Synapse
t_d_th_cor=5;
syn_func_th = const*(t_vec-t_d_th_cor).*(exp(-(t_vec-t_d_th_cor)/tau)).*((t_vec>=t_d_th_cor)&(t_vec<=t_a));

%STN-GPe Synapse
t_d_stn_gpe=2;
taudstngpea=2.5;
taurstngpea=0.4;
taudstngpen=67;
taurstngpen=2;
tpeakstngpea = t_d_stn_gpe + (((taudstngpea*taurstngpea)/(taudstngpea-taurstngpea))*log(taudstngpea/taurstngpea)); 
fstngpea = 1/(exp(-(tpeakstngpea-t_d_stn_gpe)/taudstngpea)-exp(-(tpeakstngpea-t_d_stn_gpe)/taurstngpea));
syn_func_stn_gpea = gpeak*fstngpea.*(exp(-(t_vec-t_d_stn_gpe)/taudstngpea)-exp(-(t_vec-t_d_stn_gpe)/taurstngpea)).*((t_vec>=t_d_stn_gpe)&(t_vec<=t_a));
tpeakstngpen = t_d_stn_gpe + (((taudstngpen*taurstngpen)/(taudstngpen-taurstngpen))*log(taudstngpen/taurstngpen)); 
fstngpen = 1/(exp(-(tpeakstngpen-t_d_stn_gpe)/taudstngpen)-exp(-(tpeakstngpen-t_d_stn_gpe)/taurstngpen));
syn_func_stn_gpen = gpeak*fstngpen.*(exp(-(t_vec-t_d_stn_gpe)/taudstngpen)-exp(-(t_vec-t_d_stn_gpe)/taurstngpen)).*((t_vec>=t_d_stn_gpe)&(t_vec<=t_a));

%STN-GPi Synapse
t_d_stn_gpi=1.5;
syn_func_stn_gpi = const*(t_vec-t_d_stn_gpi).*(exp(-(t_vec-t_d_stn_gpi)/tau)).*((t_vec>=t_d_stn_gpi)&(t_vec<=t_a));

%GPe-STN Synapse
t_d_gpe_stn=4;
taudg=7.7;
taurg=0.4;
tpeakg = t_d_gpe_stn + (((taudg*taurg)/(taudg-taurg))*log(taudg/taurg)); 
fg = 1/(exp(-(tpeakg-t_d_gpe_stn)/taudg)-exp(-(tpeakg-t_d_gpe_stn)/taurg));
syn_func_gpe_stn = gpeak1*fg.*(exp(-(t_vec-t_d_gpe_stn)/taudg)-exp(-(t_vec-t_d_gpe_stn)/taurg)).*((t_vec>=t_d_gpe_stn)&(t_vec<=t_a));

%GPe-GPi Synapse
t_d_gpe_gpi=3;
syn_func_gpe_gpi = const1*(t_vec-t_d_gpe_gpi).*(exp(-(t_vec-t_d_gpe_gpi)/tau)).*((t_vec>=t_d_gpe_gpi)&(t_vec<=t_a));

%GPe-GPe Synapse
t_d_gpe_gpe=1;
syn_func_gpe_gpe = const1*(t_vec-t_d_gpe_gpe).*(exp(-(t_vec-t_d_gpe_gpe)/tau)).*((t_vec>=t_d_gpe_gpe)&(t_vec<=t_a));

%GPi-TH Synapse
t_d_gpi_th=5;
syn_func_gpi_th = const1*(t_vec-t_d_gpi_th).*(exp(-(t_vec-t_d_gpi_th)/tau)).*((t_vec>=t_d_gpi_th)&(t_vec<=t_a));

%Indirect Str-GPe Synapse
t_d_d2_gpe=5;
syn_func_str_indr = const2*(t_vec- t_d_d2_gpe).*(exp(-(t_vec-t_d_d2_gpe)/tau)).*((t_vec>=t_d_d2_gpe)&(t_vec<=t_a));

%direct Str-GPi Synapse
t_d_d1_gpi=4;
syn_func_str_dr = const2*(t_vec- t_d_d1_gpi).*(exp(-(t_vec-t_d_d1_gpi)/tau)).*((t_vec>=t_d_d1_gpi)&(t_vec<=t_a));

%Cortex-Indirect Str Synapse
t_d_cor_d2=5.1;
syn_func_cor_d2 = const*(t_vec-t_d_cor_d2).*(exp(-(t_vec-t_d_cor_d2)/tau)).*((t_vec>=t_d_cor_d2)&(t_vec<=t_a));
 
%Cortex-STN Synapse
t_d_cor_stn=5.9;
taudn=90;
taurn=2;
tauda=2.49;
taura=0.5;
tpeaka = t_d_cor_stn + (((tauda*taura)/(tauda-taura))*log(tauda/taura)); 
fa = 1/(exp(-(tpeaka-t_d_cor_stn)/tauda)-exp(-(tpeaka-t_d_cor_stn)/taura));
syn_func_cor_stn_a = gpeak*fa.*(exp(-(t_vec-t_d_cor_stn)/tauda)-exp(-(t_vec-t_d_cor_stn)/taura)).*((t_vec>=t_d_cor_stn)&(t_vec<=t_a));
tpeakn = t_d_cor_stn + (((taudn*taurn)/(taudn-taurn))*log(taudn/taurn)); 
fn = 1/(exp(-(tpeakn-t_d_cor_stn)/taudn)-exp(-(tpeakn-t_d_cor_stn)/taurn));
syn_func_cor_stn_n = gpeak*fn.*(exp(-(t_vec-t_d_cor_stn)/taudn)-exp(-(t_vec-t_d_cor_stn)/taurn)).*((t_vec>=t_d_cor_stn)&(t_vec<=t_a));

t_list_th(1:n) = struct('times',[]);
t_list_cor(1:n) = struct('times',[]);
t_list_str_indr(1:n) = struct('times',[]);
t_list_str_dr(1:n) = struct('times',[]);
t_list_stn(1:n) = struct('times',[]);
t_list_gpe(1:n) = struct('times',[]);
t_list_gpi(1:n) = struct('times',[]);


all=randsample(n,n);
bll=randsample(n,n);
cll=randsample(n,n);
dll=randsample(n,n);
ell=randsample(n,n);
fll=randsample(n,n);
gll=randsample(n,n);
hll=randsample(n,n);
ill=randsample(n,n);
jll=randsample(n,n);
kll=randsample(n,n);
lll=randsample(n,n);
mll=randsample(n,n);
nll=randsample(n,n);
oll=randsample(n,n);

gcorsna=0.3*rand(n,1);
gcorsnn=0.003*rand(n,1);
gcordrstr=(0.07-0.044*pd)+0.001*rand(n,1);
ggege=1*rand(n,1);


gsngen=zeros(n,1);
gsngen(randperm(10,2)')=0.002*rand(2,1);
gsngea=zeros(n,1);
gsngea(randperm(10,2)')=0.3*rand(2,1);
gsngi=zeros(n,1);
gsngi(randperm(10,5)')=0.15;
ggith=0.112;
ggesn=0.5;
gstrgpe=0.5;
gstrgpi=0.5;
ggigi=0.5;
gm=1;
ggaba=0.1;
gcorindrstr=0.07;
gie=0.2;
gthcor=0.15;
gei=0.1;




       
    for i=2:length(t)  
        
        V1=vth(:,i-1);   
        V2=vsn(:,i-1);     
        V3=vge(:,i-1);    
        V4=vgi(:,i-1);
        V5=vstr_indr(:,i-1);
        V6=vstr_dr(:,i-1);
        V7=ve(:,i-1);
        V8=vi(:,i-1);

    % Synapse parameters 
    
    S21a(2:n)=S2a(1:n-1);
    S21a(1)=S2a(n);
    
    S21an(2:n)=S2an(1:n-1);
    S21an(1)=S2an(n);
    
    S21b(2:n)=S2b(1:n-1);
    S21b(1)=S2b(n);

    S31a(1:n-1)=S3a(2:n);
    S31a(n)=S3a(1);
    
    S31b(1:n-1)=S3b(2:n);
    S31b(n)=S3b(1);
    
    S31c(1:n-1)=S3c(2:n);
    S31c(n)=S3c(1);
    
    S32c(3:n)=S3c(1:n-2);
    S32c(1:2)=S3c(n-1:n);
    
    S32b(3:n)=S3b(1:n-2);
    S32b(1:2)=S3b(n-1:n);
    
    S11cr=S1c(all);
    S12cr=S1c(bll);
    S13cr=S1c(cll);
    S14cr=S1c(dll);

    S11br=S1b(ell);
    S12br=S1b(fll);
    S13br=S1b(gll);
    S14br=S1b(hll);

    S11ar=S1a(ill);
    S12ar=S1a(jll);
    S13ar=S1a(kll);
    S14ar=S1a(lll);

    S81r=S8(mll);
    S82r=S8(nll);
    S83r=S8(oll);

    S51(1:n-1)=S5(2:n);
    S51(n)=S5(1);
    S52(1:n-2)=S5(3:n);
    S52(n-1:n)=S5(1:2);
    S53(1:n-3)=S5(4:n);
    S53(n-2:n)=S5(1:3);
    S54(1:n-4)=S5(5:n);
    S54(n-3:n)=S5(1:4);
    S55(1:n-5)=S5(6:n);
    S55(n-4:n)=S5(1:5);
    S56(1:n-6)=S5(7:n);
    S56(n-5:n)=S5(1:6);
    S57(1:n-7)=S5(8:n);
    S57(n-6:n)=S5(1:7);
    S58(1:n-8)=S5(9:n);
    S58(n-7:n)=S5(1:8);
    S59(1:n-9)=S5(10:n);
    S59(n-8:n)=S5(1:9);

    S61b(1:n-1)=S6b(2:n);
    S61b(n)=S6b(1);
    
    S61bn(1:n-1)=S6bn(2:n);
    S61bn(n)=S6bn(1);
    
    S91(1:n-1)=S9(2:n);
    S91(n)=S9(1);
    S92(1:n-2)=S9(3:n);
    S92(n-1:n)=S9(1:2);
    S93(1:n-3)=S9(4:n);
    S93(n-2:n)=S9(1:3);
    S94(1:n-4)=S9(5:n);
    S94(n-3:n)=S9(1:4);
    S95(1:n-5)=S9(6:n);
    S95(n-4:n)=S9(1:5);
    S96(1:n-6)=S9(7:n);
    S96(n-5:n)=S9(1:6);
    S97(1:n-7)=S9(8:n);
    S97(n-6:n)=S9(1:7);
    S98(1:n-8)=S9(9:n);
    S98(n-7:n)=S9(1:8);
    S99(1:n-9)=S9(10:n);
    S99(n-8:n)=S9(1:9);
    
    m1=th_minf(V1);
    m3=gpe_minf(V3);m4=gpe_minf(V4);
    n3=gpe_ninf(V3);n4=gpe_ninf(V4);
    h1=th_hinf(V1);
    h3=gpe_hinf(V3);h4=gpe_hinf(V4);
    p1=th_pinf(V1);
    a3=gpe_ainf(V3);a4=gpe_ainf(V4);
    s3=gpe_sinf(V3);s4=gpe_sinf(V4);
    r1=th_rinf(V1);
    r3=gpe_rinf(V3);r4=gpe_rinf(V4);

    tn3=gpe_taun(V3);tn4=gpe_taun(V4);
    th1=th_tauh(V1);
    th3=gpe_tauh(V3);th4=gpe_tauh(V4);
    tr1=th_taur(V1);tr3=30;tr4=30;
    
    n2=stn_ninf(V2);
    m2=stn_minf(V2);
    h2=stn_hinf(V2);
    a2=stn_ainf(V2);
    b2=stn_binf(V2);
    c2=stn_cinf(V2);
    d2=stn_d2inf(V2);
    d1=stn_d1inf(V2);
    p2=stn_pinf(V2);
    q2=stn_qinf(V2);
    r2=stn_rinf(V2);
 
    td2=130;
    tr2=2;
    tn2=stn_taun(V2);
    tm2=stn_taum(V2);
    th2=stn_tauh(V2);
    ta2=stn_taua(V2);
    tb2=stn_taub(V2);
    tc2=stn_tauc(V2);
    td1=stn_taud1(V2);
    tp2=stn_taup(V2);
    tq2=stn_tauq(V2);

    Ecasn=con*log(Cao./CAsn2);
   
   
    %thalamic cell currents
    Il1=gl(1)*(V1-El(1));
    Ina1=gna(1)*(m1.^3).*H1.*(V1-Ena(1));
    Ik1=gk(1)*((0.75*(1-H1)).^4).*(V1-Ek(1));
    It1=gt(1)*(p1.^2).*R1.*(V1-Et);
    Igith=ggith*(V1-Esyn(6)).*(S4); 
    Iappth=1.2;
 

    %STN cell currents
    Ina2=gna(2)*(M2.^3).*H2.*(V2-Ena(2));
    Ik2=gk(2)*(N2.^4).*(V2-Ek(2));
    Ia2=ga*(A2.^2).*(B2).*(V2-Ek(2));
    IL2=gL*(C2.^2).*(D1).*(D2).*(V2-Ecasn);
    It2=(gt(2)*(P2.^2).*(Q2).*(V2-Ecasn));
    Icak2=gcak*(R2.^2).*(V2-Ek(2));
    Il2=gl(2)*(V2-El(2));
    Igesn=(ggesn*((V2-Esyn(1)).*(S3a+S31a))); 
    Icorsnampa=gcorsna.*(V2-Esyn(2)).*(S6b+S61b);
    Icorsnnmda=gcorsnn.*(V2-Esyn(2)).*(S6bn+S61bn);

    %GPe cell currents
    Il3=gl(3)*(V3-El(3));
    Ik3=gk(3)*(N3.^4).*(V3-Ek(3));
    Ina3=gna(3)*(m3.^3).*H3.*(V3-Ena(3));
    It3=gt(3)*(a3.^3).*R3.*(V3-Eca(3));
    Ica3=gca(3)*(s3.^2).*(V3-Eca(3));
    Iahp3=gahp(3)*(V3-Ek(3)).*(CA3./(CA3+k1(3)));
    Isngeampa=(gsngea).*((V3-Esyn(2)).*(S2a+S21a)); 
    Isngenmda=(gsngen).*((V3-Esyn(2)).*(S2an+S21an)); 
    Igege=(0.25*(pd*3+1))*(ggege).*((V3-Esyn(3)).*(S31c+S32c)); 
    Istrgpe=gstrgpe*(V3-Esyn(6)).*(S5+S51+S52+S53+S54+S55+S56+S57+S58+S59);
    Iappgpe=3-2*corstim*~pd; %Modulation only during cortical stim to maintain mean firing rate

    %GPi cell currents
    Il4=gl(3)*(V4-El(3));
    Ik4=gk(3)*(N4.^4).*(V4-Ek(3));
    Ina4=gna(3)*(m4.^3).*H4.*(V4-Ena(3));
    It4=gt(3)*(a4.^3).*R4.*(V4-Eca(3));
    Ica4=gca(3)*(s4.^2).*(V4-Eca(3));
    Iahp4=gahp(3)*(V4-Ek(3)).*(CA4./(CA4+k1(3)));
    Isngi=(gsngi).*((V4-Esyn(4)).*(S2b+S21b));
    Igigi=ggigi*((V4-Esyn(5)).*(S31b+S32b)); 
    Istrgpi=gstrgpi*(V4-Esyn(6)).*(S9+S91+S92+S93+S94+S95+S96+S97+S98+S99);
    Iappgpi=3;

    %Striatum D2 cell currents
    Ina5=gna(4)*(m5.^3).*h5.*(V5-Ena(4));
    Ik5=gk(4)*(n5.^4).*(V5-Ek(4));
    Il5=gl(4)*(V5-El(4));
    Im5=(2.6-1.1*pd)*gm*p5.*(V5-Em);
    Igaba5=(ggaba/4)*(V5-Esyn(7)).*(S11cr+S12cr+S13cr+S14cr);
    Icorstr5=gcorindrstr*(V5-Esyn(2)).*(S6a);
    
    %Striatum D1 cell currents
    Ina6=gna(4)*(m6.^3).*h6.*(V6-Ena(4));
    Ik6=gk(4)*(n6.^4).*(V6-Ek(4));
    Il6=gl(4)*(V6-El(4));
    Im6=(2.6-1.1*pd)*gm*p6.*(V6-Em);
    Igaba6=(ggaba/3)*(V6-Esyn(7)).*(S81r+S82r+S83r);
    Icorstr6=gcordrstr.*(V6-Esyn(2)).*(S6a);
    
    %Excitatory Neuron Currents
    Iie=gie*(V7-Esyn(1)).*(S11br+S12br+S13br+S14br);
    Ithcor=gthcor*(V7-Esyn(2)).*(S7);

    
    %Inhibitory Neuron Currents
    Iei=gei*(V8-Esyn(2)).*(S11ar+S12ar+S13ar+S14ar);

    %Differential Equations for cells
    %thalamic
    vth(:,i)= V1+dt*(1/Cm*(-Il1-Ik1-Ina1-It1-Igith+Iappth));
    H1=H1+dt*((h1-H1)./th1);
    R1=R1+dt*((r1-R1)./tr1);

for j=1:n

if (vth(j,i-1)<-10 && vth(j,i)>-10) % check for input spike
     t_list_th(j).times = [t_list_th(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S7(j) = sum(syn_func_th(t_list_th(j).times));

   % Update spike times
   if t_list_th(j).times
     t_list_th(j).times = t_list_th(j).times + 1;
     if (t_list_th(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_th(j).times = t_list_th(j).times((2:max(size(t_list_th(j).times))));
     end
   end
end
    
    %STN

    vsn(:,i)=V2+dt*(1/Cm*(-Ina2-Ik2-Ia2-IL2-It2-Icak2-Il2-Igesn-Icorsnampa-Icorsnnmda+Idbs(i))); %STN-DBS
    N2=N2+dt*((n2-N2)./tn2); 
    H2=H2+dt*((h2-H2)./th2);
    M2=M2+dt*((m2-M2)./tm2); 
    A2=A2+dt*((a2-A2)./ta2);
    B2=B2+dt*((b2-B2)./tb2); 
    C2=C2+dt*((c2-C2)./tc2);
    D2=D2+dt*((d2-D2)./td2); 
    D1=D1+dt*((d1-D1)./td1);
    P2=P2+dt*((p2-P2)./tp2); 
    Q2=Q2+dt*((q2-Q2)./tq2);
    R2=R2+dt*((r2-R2)./tr2); 
    
    CAsn2=CAsn2+dt*((-alp*(IL2+It2))-(Kca*CAsn2));
for j=1:n

if (vsn(j,i-1)<-10 && vsn(j,i)>-10) % check for input spike
     t_list_stn(j).times = [t_list_stn(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S2a(j) = sum(syn_func_stn_gpea(t_list_stn(j).times));
   S2an(j) = sum(syn_func_stn_gpen(t_list_stn(j).times));

   S2b(j) = sum(syn_func_stn_gpi(t_list_stn(j).times));

   % Update spike times
   if t_list_stn(j).times
     t_list_stn(j).times = t_list_stn(j).times + 1;
     if (t_list_stn(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_stn(j).times = t_list_stn(j).times((2:max(size(t_list_stn(j).times))));
     end
   end
end
    
    %GPe
    vge(:,i)=V3+dt*(1/Cm*(-Il3-Ik3-Ina3-It3-Ica3-Iahp3-Isngeampa-Isngenmda-Igege-Istrgpe+Iappgpe));
    N3=N3+dt*(0.1*(n3-N3)./tn3);
    H3=H3+dt*(0.05*(h3-H3)./th3);
    R3=R3+dt*(1*(r3-R3)./tr3);
    CA3=CA3+dt*(1*10^-4*(-Ica3-It3-kca(3)*CA3));
for j=1:n

if (vge(j,i-1)<-10 && vge(j,i)>-10) % check for input spike
     t_list_gpe(j).times = [t_list_gpe(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S3a(j) = sum(syn_func_gpe_stn(t_list_gpe(j).times));
   S3b(j) = sum(syn_func_gpe_gpi(t_list_gpe(j).times));
   S3c(j) = sum(syn_func_gpe_gpe(t_list_gpe(j).times));


   % Update spike times
   if t_list_gpe(j).times
     t_list_gpe(j).times = t_list_gpe(j).times + 1;
     if (t_list_gpe(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_gpe(j).times = t_list_gpe(j).times((2:max(size(t_list_gpe(j).times))));
     end
   end
end
    
    %GPi
    vgi(:,i)=V4+dt*(1/Cm*(-Il4-Ik4-Ina4-It4-Ica4-Iahp4-Isngi-Igigi-Istrgpi+Iappgpi));
    N4=N4+dt*(0.1*(n4-N4)./tn4);
    H4=H4+dt*(0.05*(h4-H4)./th4);
    R4=R4+dt*(1*(r4-R4)./tr4);
    CA4=CA4+dt*(1*10^-4*(-Ica4-It4-kca(3)*CA4));

for j=1:n

if (vgi(j,i-1)<-10 && vgi(j,i)>-10) % check for input spike
     t_list_gpi(j).times = [t_list_gpi(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S4(j) = sum(syn_func_gpi_th(t_list_gpi(j).times));


   % Update spike times
   if t_list_gpi(j).times
     t_list_gpi(j).times = t_list_gpi(j).times + 1;
     if (t_list_gpi(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_gpi(j).times = t_list_gpi(j).times((2:max(size(t_list_gpi(j).times))));
     end
   end
end
    
    %Striatum D2
 vstr_indr(:,i)=V5+(dt/Cm)*(-Ina5-Ik5-Il5-Im5-Igaba5-Icorstr5);
 m5=m5+dt*(alpham(V5).*(1-m5)-betam(V5).*m5);
 h5=h5+dt*(alphah(V5).*(1-h5)-betah(V5).*h5);
 n5=n5+dt*(alphan(V5).*(1-n5)-betan(V5).*n5);
 p5=p5+dt*(alphap(V5).*(1-p5)-betap(V5).*p5);
 S1c=S1c+dt*((Ggaba(V5).*(1-S1c))-(S1c/tau_i));

for j=1:n

if (vstr_indr(j,i-1)<-10 && vstr_indr(j,i)>-10) % check for input spike
     t_list_str_indr(j).times = [t_list_str_indr(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S5(j) = sum(syn_func_str_indr(t_list_str_indr(j).times));

   % Update spike times
   if t_list_str_indr(j).times
     t_list_str_indr(j).times = t_list_str_indr(j).times + 1;
     if (t_list_str_indr(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_str_indr(j).times = t_list_str_indr(j).times((2:max(size(t_list_str_indr(j).times))));
     end
   end
end

% %Striatum D1 type
 vstr_dr(:,i)=V6+(dt/Cm)*(-Ina6-Ik6-Il6-Im6-Igaba6-Icorstr6);
 m6=m6+dt*(alpham(V6).*(1-m6)-betam(V6).*m6);
 h6=h6+dt*(alphah(V6).*(1-h6)-betah(V6).*h6);
 n6=n6+dt*(alphan(V6).*(1-n6)-betan(V6).*n6);
 p6=p6+dt*(alphap(V6).*(1-p6)-betap(V6).*p6);
 S8=S8+dt*((Ggaba(V6).*(1-S8))-(S8/tau_i));

 
 for j=1:n

if (vstr_dr(j,i-1)<-10 && vstr_dr(j,i)>-10) % check for input spike
     t_list_str_dr(j).times = [t_list_str_dr(j).times; 1];
end   
   % Calculate synaptic current due to current and past input spikes
   S9(j) = sum(syn_func_str_dr(t_list_str_dr(j).times));

   % Update spike times
   if t_list_str_dr(j).times
     t_list_str_dr(j).times = t_list_str_dr(j).times + 1;
     if (t_list_str_dr(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_str_dr(j).times = t_list_str_dr(j).times((2:max(size(t_list_str_dr(j).times))));
     end
   end
 end

%Excitatory Neuron
    ve(:,i)=V7+dt*((0.04*(V7.^2))+(5*V7)+140-ue(:,i-1)-Iie-Ithcor+Iappco(i));
    ue(:,i)=ue(:,i-1)+dt*(ae*((be*V7)-ue(:,i-1)));
    
   for j=1:n
        if ve(j,i-1)>=30
        ve(j,i)=ce;
        ue(j,i)=ue(j,i-1)+de;
        
 t_list_cor(j).times = [t_list_cor(j).times; 1];
        end
   
   % Calculate synaptic current due to current and past input spikes
   S6a(j) = sum(syn_func_cor_d2(t_list_cor(j).times));
   S6b(j) = sum(syn_func_cor_stn_a(t_list_cor(j).times));
   S6bn(j) = sum(syn_func_cor_stn_n(t_list_cor(j).times));

   % Update spike times
   if t_list_cor(j).times
     t_list_cor(j).times = t_list_cor(j).times + 1;
     if (t_list_cor(j).times(1) == t_a/dt)  % Reached max duration of syn conductance
       t_list_cor(j).times = t_list_cor(j).times((2:max(size(t_list_cor(j).times))));
     end
   end
   
   end        
    
    ace=find(ve(:,i-1)<-10 & ve(:,i)>-10);
    uce=zeros(n,1); uce(ace)=gpeak/(tau*exp(-1))/dt;
    S1a=S1a+dt*Z1a; 
    z1adot=uce-2/tau*Z1a-1/(tau^2)*S1a;
    Z1a=Z1a+dt*z1adot;
    
    %Inhibitory InterNeuron
    vi(:,i)=V8+dt*((0.04*(V8.^2))+(5*V8)+140-ui(:,i-1)-Iei+Iappco(i));
    ui(:,i)=ui(:,i-1)+dt*(ai*((bi*V8)-ui(:,i-1)));
    
   for j=1:n
        if vi(j,i-1)>=30
        vi(j,i)=ci;
        ui(j,i)=ui(j,i-1)+di;
        end
   end
        
    
    aci=find(vi(:,i-1)<-10 & vi(:,i)>-10);
    uci=zeros(n,1); uci(aci)=gpeak/(tau*exp(-1))/dt;
    S1b=S1b+dt*Z1b; 
    z1bdot=uci-2/tau*Z1b-1/(tau^2)*S1b;
    Z1b=Z1b+dt*z1bdot;


    end

    
    [TH_APs]  = find_spike_times(vth,t,n);
    [STN_APs] = find_spike_times(vsn,t,n);
    [GPe_APs] = find_spike_times(vge,t,n);
    [GPi_APs] = find_spike_times(vgi,t,n);
    [Striat_APs_indr]=find_spike_times(vstr_indr,t,n);
    [Striat_APs_dr]=find_spike_times(vstr_dr,t,n);
    [Cor_APs] = find_spike_times([ve;vi],t,2*n);
end

function [data] = find_spike_times(v,t,nn)


    data(1:nn) = struct('times',[]);
    t = t./1000;    % unit: second
    for k = 1:nn
        data(k).times = t(diff(v(k,:)>-20)==1)';
    end

end

function [ainf] = gpe_ainf(V)
    ainf=1./(1+exp(-(V+57)./2));
end

function [hinf] = gpe_hinf(V)
    hinf=1./(1+exp((V+58)./12));
end

function [minf] = gpe_minf(V)
    minf=1./(1+exp(-(V+37)./10));
end

function [ninf] = gpe_ninf(V)
    ninf=1./(1+exp(-(V+50)./14));
end

function [rinf] = gpe_rinf(V)
    rinf=1./(1+exp((V+70)./2));
end

function [sinf] = gpe_sinf(V)
    sinf=1./(1+exp(-(V+35)./2));
end

function [tau] = gpe_tauh(V)
    tau=0.05+0.27./(1+exp(-(V+40)./-12));
end

function [tau] = gpe_taun(V)
    tau=0.05+0.27./(1+exp(-(V+40)./-12));
end

function [hinf] = th_hinf(V)
    hinf=1./(1+exp((V+41)./4));
end

function [minf] = th_minf(V)
    minf=1./(1+exp(-(V+37)./7));
end

function [pinf] = th_pinf(V)
    pinf=1./(1+exp(-(V+60)./6.2));
end

function [rinf] = th_rinf(V)
    rinf=1./(1+exp((V+84)./4));
end

function [tau] = th_tauh(V)
    tau=1./(ah(V)+bh(V));
end

function [a] = ah(V)
    a=0.128*exp(-(V+46)./18); % part of th_tauh fxn
end

function [b] = bh(V)
    b=4./(1+exp(-(V+23)./5)); % part of th_tauh fxn
end

function [tau] = th_taur(V)
    tau=0.15*(28+exp(-(V+25)./10.5));
end

function [ah] = alphah(V)
ah=0.128*exp((-50-V)/18);
end

function [am] = alpham(V)
am=(0.32*(54+V))./(1-exp((-54-V)/4));
end

function [an] = alphan(V)
an=(0.032*(52+V))./(1-exp((-52-V)./5));
end

function [ap] = alphap(V)
ap=(3.209*10^-4*(30+V))./(1-exp((-30-V)./9));
end

function [bh] = betah(V)
bh=4./(1+exp((-27-V)/5));
end

function [bn] = betan(V)
bn=0.5*exp((-57-V)./40);
end

function [bm] = betam(V)
bm=0.28*(V+27)./((exp((27+V)/5))-1);
end

function [bp] = betap(V)
bp=(-3.209*10^-4*(30+V))./(1-exp((30+V)./9));
end

function [gb] = Ggaba(V)
gb=2*(1+tanh(V/4));
end

function [ainf] = stn_ainf(V)
    ainf=1./(1+exp(-(V+45)./14.7));
end

function [binf] = stn_binf(V)
    binf=1./(1+exp((V+90)./7.5));
end

function [cinf] = stn_cinf(V)
    cinf=1./(1+exp(-(V+30.6)./5));
end

function [d1inf] = stn_d1inf(V)
    d1inf=1./(1+exp((V+60)./7.5));
end

function [d2inf] = stn_d2inf(V)
    d2inf=1./(1+exp((V-0.1)./0.02));
end

function [hinf] = stn_hinf(V)
    hinf=1./(1+exp((V+45.5)./6.4));
end

function [minf] = stn_minf(V)
    minf=1./(1+exp(-(V+40)./8));
end

function [ninf] = stn_ninf(V)
    ninf=1./(1+exp(-(V+41)./14));
end

function [pinf] = stn_pinf(V)
    pinf=1./(1+exp(-(V+56)./6.7));
end

function [qinf] = stn_qinf(V)
    qinf=1./(1+exp((V+85)./5.8));
end

function [rinf] = stn_rinf(V)
    rinf=1./(1+exp(-(V-0.17)./0.08));
end

function [tau] = stn_taua(V)
    tau=1+1./(1+exp(-(V+40)./-0.5));
end

function [tau] = stn_taub(V)
    tau=200./(exp(-(V+60)./-30)+exp(-(V+40)./10));
end

function [tau] = stn_tauc(V)
    tau=45+10./(exp(-(V+27)./-20)+exp(-(V+50)./15));
end

function [tau] = stn_taud1(V)
    tau=400+500./(exp(-(V+40)./-15)+exp(-(V+20)./20));
end

function [tau] = stn_tauh(V)
    tau=24.5./(exp(-(V+50)./-15)+exp(-(V+50)./16));
end

function [tau] = stn_taum(V)
    tau=0.2+3./(1+exp(-(V+53)./-0.7));
end

function [tau] = stn_taun(V)
    tau=11./(exp(-(V+40)./-40)+exp(-(V+40)./50));
end

function [tau] = stn_taup(V)
    tau=5+0.33./(exp(-(V+27)./-10)+exp(-(V+102)./15));
end

function [tau] = stn_tauq(V)
    tau=400./(exp(-(V+50)./-15)+exp(-(V+50)./16));
end

function [area S f] = make_Spectrum(raw,params)

% Compute Multitaper Spectrum
[S,f] = mtspectrumpt(raw,params);
beta = S(f>7 & f<35);
betaf = f(f>7 & f<35);
area = trapz(betaf,beta);

end

function [S,f,R,Serr]=mtspectrumpt(data,params,fscorr,t)

if nargin < 1; error('Need data'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear params
data=change_row_to_column(data);
if nargout > 3 && err(1)==0; error('cannot compute error bars with err(1)=0; change params and run again'); end;
if nargin < 3 || isempty(fscorr); fscorr=0;end;
if nargin < 4 || isempty(t);
   [mintime,maxtime]=minmaxsptimes(data);
   dt=1/Fs; % sampling time
   t=mintime-dt:dt:maxtime+dt; % time grid for prolates
end;
N=length(t); % number of points in grid for dpss
nfft=max(2^(nextpow2(N)+pad),N); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers=dpsschk(tapers,N,Fs); % check tapers
[J,Msp,Nsp]=mtfftpt(data,tapers,nfft,t,f,findx); % mt fft for point process times
S=squeeze(mean(conj(J).*J,2));
if trialave; S=squeeze(mean(S,2));Msp=mean(Msp);end;
R=Msp*Fs;
if nargout==4;
   if fscorr==1;
      Serr=specerr(S,J,err,trialave,Nsp);
   else
      Serr=specerr(S,J,err,trialave);
   end;
end;
end

function [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)

if ~isfield(params,'tapers') || isempty(params.tapers);  %If the tapers don't exist
     display('tapers unspecified, defaulting to params.tapers=[3 5]');
     params.tapers=[3 5];
end;
if ~isempty(params) && length(params.tapers)==3 
    % Compute timebandwidth product
    TW = params.tapers(2)*params.tapers(1);
    % Compute number of tapers
    K  = floor(2*TW - params.tapers(3));
    params.tapers = [TW  K];
end

if ~isfield(params,'pad') || isempty(params.pad);
    params.pad=0;
end;
if ~isfield(params,'Fs') || isempty(params.Fs);
    params.Fs=1;
end;
if ~isfield(params,'fpass') || isempty(params.fpass);
    params.fpass=[0 params.Fs/2];
end;
if ~isfield(params,'err') || isempty(params.err);
    params.err=0;
end;
if ~isfield(params,'trialave') || isempty(params.trialave);
    params.trialave=0;
end;

tapers=params.tapers;
pad=params.pad;
Fs=params.Fs;
fpass=params.fpass;
err=params.err;
trialave=params.trialave;
end

function data=change_row_to_column(data)

dtmp=[];
if isstruct(data);
   C=length(data);
   if C==1;
      fnames=fieldnames(data);
      eval(['dtmp=data.' fnames{1} ';'])
      data=dtmp(:);
   end
else
  [N,C]=size(data);
  if N==1 || C==1;
    data=data(:);
  end;
end;
end

function [mintime, maxtime]=minmaxsptimes(data)

dtmp='';
if isstruct(data)
   data=reshape(data,numel(data),1);
   C=size(data,1);
   fnames=fieldnames(data);
   mintime=zeros(1,C); maxtime=zeros(1,C);
   for ch=1:C
     eval(['dtmp=data(ch).' fnames{1} ';'])
     if ~isempty(dtmp)
        maxtime(ch)=max(dtmp);
        mintime(ch)=min(dtmp);
     else
        mintime(ch)=NaN;
        maxtime(ch)=NaN;
     end
   end;
   maxtime=max(maxtime); % maximum time
   mintime=min(mintime); % minimum time
else
     dtmp=data;
     if ~isempty(dtmp)
        maxtime=max(dtmp);
        mintime=min(dtmp);
     else
        mintime=NaN;
        maxtime=NaN;
     end
end
if mintime < 0 
   error('Minimum spike time is negative'); 
end
end

function [f,findx]=getfgrid(Fs,nfft,fpass)

if nargin < 3; error('Need all arguments'); end;
df=Fs/nfft;
f=0:df:Fs; % all possible frequencies
f=f(1:nfft);
if length(fpass)~=1;
   findx=find(f>=fpass(1) & f<=fpass(end));
else
   [fmin,findx]=min(abs(f-fpass));
   clear fmin
end;
f=f(findx);
end

function [tapers,eigs]=dpsschk(tapers,N,Fs)

if nargin < 3; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 && sz(2)==2;
    [tapers,eigs]=dpss(N,tapers(1),tapers(2));
    tapers = tapers*sqrt(Fs);
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;
end

function [J,Msp,Nsp]=mtfftpt(data,tapers,nfft,t,f,findx)

if nargin < 6; error('Need all input arguments'); end;
if isstruct(data); C=length(data); else C=1; end% number of channels
K=size(tapers,2); % number of tapers
nfreq=length(f); % number of frequencies
if nfreq~=length(findx); error('frequency information (last two arguments) inconsistent'); end;
H=fft(tapers,nfft,1);  % fft of tapers
H=H(findx,:); % restrict fft of tapers to required frequencies
w=2*pi*f; % angular frequencies at which ft is to be evaluated
Nsp=zeros(1,C); Msp=zeros(1,C);
for ch=1:C;
  if isstruct(data);
     fnames=fieldnames(data);
     eval(['dtmp=data(ch).' fnames{1} ';'])
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  else
     dtmp=data;
     indx=find(dtmp>=min(t)&dtmp<=max(t));
     if ~isempty(indx); dtmp=dtmp(indx);
     end;
  end;
  Nsp(ch)=length(dtmp);
  Msp(ch)=Nsp(ch)/length(t);
  if Msp(ch)~=0;
      data_proj=interp1(t',tapers,dtmp);
      exponential=exp(-1i*w'*(dtmp-t(1))');
      J(:,:,ch)=exponential*data_proj-H*Msp(ch);
  else
      J(1:nfreq,1:K,ch)=0;
  end;
end;
end

function Serr=specerr(S,J,err,trialave,numsp)
 
if nargin < 4; error('Need at least 4 input arguments'); end;
if err(1)==0; error('Need err=[1 p] or [2 p] for error bar calculation. Make sure you are not asking for the output of Serr'); end;
[nf,K,C]=size(J);
errchk=err(1);
p=err(2);
pp=1-p/2;
qq=1-pp;

if trialave
   dim=K*C;
   C=1;
   dof=2*dim;
   if nargin==5; dof = fix(1/(1/dof + 1/(2*sum(numsp)))); end
   J=reshape(J,nf,dim);
else
   dim=K;
   dof=2*dim*ones(1,C);
   for ch=1:C;
     if nargin==5; dof(ch) = fix(1/(1/dof + 1/(2*numsp(ch)))); end 
   end;
end;
Serr=zeros(2,nf,C);
if errchk==1;
   Qp=chi2inv(pp,dof);
   Qq=chi2inv(qq,dof);
   Serr(1,:,:)=dof(ones(nf,1),:).*S./Qp(ones(nf,1),:);
   Serr(2,:,:)=dof(ones(nf,1),:).*S./Qq(ones(nf,1),:);
elseif errchk==2;
   tcrit=tinv(pp,dim-1);
   for k=1:dim;
       indices=setdiff(1:dim,k);
       Jjk=J(:,indices,:); % 1-drop projection
       eJjk=squeeze(sum(Jjk.*conj(Jjk),2));
       Sjk(k,:,:)=eJjk/(dim-1); % 1-drop spectrum
   end;
   sigma=sqrt(dim-1)*squeeze(std(log(Sjk),1,1)); if C==1; sigma=sigma'; end; 
   conf=repmat(tcrit,nf,C).*sigma;
   conf=squeeze(conf); 
   Serr(1,:,:)=S.*exp(-conf); Serr(2,:,:)=S.*exp(conf);
end;
Serr=squeeze(Serr);
end

