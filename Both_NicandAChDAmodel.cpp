#include <stdio.h>
#include <stdlib.h>
#include "mex.h" // Standard include library for the MEX-files.
 #include "matrix.h" // Standard include library for accessing MATLAB data struct.
 #include <string.h> 
#include <math.h>
#include <time.h>
 #include <iostream>

//---------- parameters for 1-compartmental model------------
double Cgaba=10;
double gL, gbarK=1/2, gbarCa=2.5/2, gbarKCa=7.8/2, gbarDR=0*10, gbarERG=0*4.8, gNa=0,gSNa=0.13; //conductance amplitudes
double gDR, gERG, gCa, gKCa, gK; //conductances
double gampa, gnmda, gbarnmda, ggaba; //synaptic conductances
double nmdasig, ampasig, allgaba; //synaptic activation functions
double k=160, buf=0.00023, zF=0.019298, tc=0.52, r=0.2; //for Ca equation
double ergtaumax=300, ergtaumin=62, Vnh=-50.4, Sn=1.8, sh=13; //for ERG
double Mg=0.5, me=0.062; //for gnmda 
double EL=-35, EK=-90, ECa=50, ENMDA=0, EAMPA=0, EGABA=-90, ENa=55; // reversal potentials
double VHK=-10, VSK=7; //for gK
double th=0.05; //for gNa
double aact=1., tades=6.1, adesrel=40, adeact=1.6, nact=7., ndeact=170.; //for AMPA 
//double nmdathresh=0.2, nmdasl=0.03, ampathresh=0.2, ampasl=0.03; // old for NMDA
double nmdathresh=9, nmdasl=1.3, ampathresh=9, ampasl=1.3; //for NMDA
double TT, dt=0.02, tgl; // time, step, global time
int N, NG=20, Neq=15+NG*4; // N of steps,  N of GABA neurons, N of equations,
double amg, bmg, minfgg, ang, bng, ahg, bhg, gnmdagg, gspikeg, vgnz; // for GABA population
double glg=0.1*10, gna=22*10, gdrg=7*10, tng=1/2.4, thg=5, tbn=0.7, as=12, bs=0.1; // for GABA population
double gmf, dg=0; //strength of Gap junction
double gampag=0, gnmdabarg=0; //GABA synaptic conductances
double inp, inpachda, inpachgaba, inpnic; // input
double vhna=-50, slna=5; // for subthresh Na current
double Vcah=-52, Sca=3; // for Ca current
double caleak=0.1; // fraction of Ca current in a leak current
double vhh=-95, slh=8, th0=625, vtauh=-112, Eh=-20, gh, gbarh; // Ih (Okamoto)
//--------------for Ach a2b2 receptor from paper---------------
double EC50_actB2=30, Hill_actB2=1.05, w_actNicB2_da=3, w_actNicB2_gaba=3;
double EC50_desB2=0.061; // half-maximum concentration of desensitization by Nic (microM)
double Hill_desB2=0.5; // Hill coefficient of desensitization
double tmin_desB2=500; // minimal des. time constant (500 msec)
double tmax_desB2=600*10^3; // maximal des. time constant (10 min)
double Ktau_desB2=0.11; // half maximum conc of desensitization time constant (microM)
double Htau_desB2=3; // Hill coefficient of des. time constant
double gachgaba, gachda, gbarachgaba, gbarachda;
double EAch=0; // reversal potential
double achtau, achinfda, achinfgaba;
double achtaudes, achinfdes;
//----------------------------------------

double *y;
double *f;
double *rngL=new double[NG];
double *glg1=new double[NG];

void system (void)
{    
double alphan=-0.0032*(y[0]+5.)/(exp(-(y[0]+5.)/10.) - 1.);
double betan=0.05*exp(-((y[0]+10.)/16.));
double alpham=-0.32*(y[0]+39)/(exp(-(y[0]+39)/4)-1);
double betam=0.28*(y[0]+4)/(exp((y[0]+4)/5)-1);
double alphah=0.2*th*exp(-((y[0]+47.)/18.));
double betah=25.*th/(1.+(exp(-(y[0]+24.)/5.)));
double ergtau=ergtaumin + ergtaumax*(1/(1 + exp((y[0]-Vnh)/Sn))-1/(1+exp((y[0]-Vnh+sh)/Sn)));
double erginf=1/(1+exp(-(y[0]-Vnh-3)/Sn));
double nsqr=y[5]*y[5];
double minf=alpham/(alpham+betam);
gDR=gbarDR*(nsqr*nsqr);
double ergsqr=y[4]*y[4];
gERG=gbarERG*(ergsqr*ergsqr);
gK=gbarK/(1. + exp(-(y[0]-VHK)/VSK));
//double nmdasig=1/(1+exp(-(y[1]-nmdathresh)/nmdasl));
//double ampasig=1/(1+exp(-(y[2]*y[6]-ampathresh)/ampasl));
double nmdasig=1/(1+exp(-(inp-nmdathresh)/nmdasl));
double ampasig=1/(1+exp(-(inp-ampathresh)/ampasl));
double alphac=((fabs(y[0]-Vcah))>0.00001)? (-0.0032*(y[0]-Vcah)/(exp(-(y[0]-Vcah)/Sca) - 1.)) : (-0.0032*0.00001/(exp(-0.00001/Sca)-1.));
double betac=0.05*exp(-(y[0]-Vcah+5)/40.);
double csinf=alphac/(alphac+betac);
double csqr=csinf*csinf;
gCa=gbarCa*(csqr*csqr);
double ksq=k*k;
double casq=y[3]*y[3];
gKCa=gbarKCa*(casq*casq)/((casq*casq) + (ksq*ksq));
gnmda=gbarnmda/(1+0.28*Mg*exp(-me*y[0]));
double na=1/(1+exp(-(y[0]-vhna)/slna)); // subthresh Na
//achsiggaba=1/(1+exp(-(y[8]*y[9]-achthreshgaba)/achslgaba));
//achsigda=1/(1+exp(-(y[8]*y[9]-achthreshda)/achslda));
//---------Ih-------------
double hinf=1/(1+exp((y[0]-vhh)/slh));
double tauh=th0*exp(0.075*(y[0]-vtauh))/(1+exp(0.083*(y[0]-vtauh)));
//double tauh=550+1860*exp(-(pow((y[0]-(-83.6))/19.4,2)));
gh=gbarh*y[11];
//--------Iach-----------
gachda = gbarachda*y[8]*(1-y[10]);
gachgaba = gbarachgaba*y[9]*(1-y[10]);
achtau=1; // 5 ms // was 1 for some reasone ( everything worked), changed back to 5
achinfdes = 1/(1+pow(EC50_desB2/inpnic,Hill_desB2));
achtaudes = tmin_desB2 + tmax_desB2 * 1/(1+pow(inpnic/Ktau_desB2,Htau_desB2));
achinfda = 1/(1+pow(EC50_actB2/(inpachda+w_actNicB2_da*inpnic),Hill_actB2));
achinfgaba = 1/(1+pow(EC50_actB2/(inpachgaba+w_actNicB2_gaba*inpnic),Hill_actB2));

allgaba=0;
for (int ig=0; ig<NG; ig++){
	allgaba+=y[15+ig*4];
}
allgaba/=NG;

f[0]= gnmda*y[1]*(ENMDA-y[0])+ gampa*y[2]*y[6]*(EAMPA-y[0])+ ggaba*allgaba*(EGABA-y[0]) + gCa*(ECa-y[0])
+ (gKCa+gK+gERG/*+gDR*/)*(EK-y[0])+ gL*(EL-y[0])+ gNa*(ENa-y[0]) + gSNa*na*(ENa-y[0]) + gachda*(EAch-y[0]) + gh*(Eh-y[0]);// + ghNa*(ENa-y[0]) + ghK*(EK-y[0]);
f[1]=nmdasig*(1-y[1])/nact-(1-nmdasig)*y[1]/ndeact; // NMDA
f[2]=ampasig*(1-y[2])/aact-(1-ampasig)*y[2]/adeact; // AMPA
f[3]= 2.*buf*((gCa+caleak*gL)*(ECa - y[0])/zF - y[3]/tc)/r; //Ca
f[4]= (erginf-y[4])/ergtau; //ERG
f[5]= alphan*(1.-y[5])-betan*y[5]; // DR
f[6]=(1-ampasig)*(1-y[6])/adesrel-ampasig*y[6]/tades; //ampa desensitization
f[7]= alphah*(1.-y[7])-betah*y[7]; // Na
f[8]= (achinfda-y[8])/(achtau); //Ach
f[9]= (achinfgaba-y[9])/achtau; //Ach
f[10]= (achinfdes-y[10])/achtaudes; //Ach
f[11]= (hinf-y[11])/tauh; //Ih

// population of GABA neurons
gmf=0;
for (int ig=0; ig<NG; ig++){gmf+=y[12+ig*4];} //for gap junctions
gmf/=NG;

for (int ig=0; ig<NG; ig++)
{
double amg=0.1*(y[12+ig*4]+30.0)/(1.0-exp(-(y[12+ig*4]+30.0)/10.0));
double bmg=4.0*exp(-(y[12+ig*4]+55.0)/18.0);
double minfgg=amg/(amg+bmg);
double ahg=0.07*exp(-(y[12+ig*4]+53.0)/20.0);
double bhg=1.0/(1.0+exp(-(y[12+ig*4]+23.0)/10.0));
double ang=0.01*(y[12+ig*4]+29.0)/(1.0-exp(-(y[12+ig*4]+29.0)/10.0));
double bng=tbn*0.125*exp(-(y[12+ig*4]+39.0)/80.0);
double gnmdagg=gnmdabarg/(1+0.28*Mg*exp(-me*y[12+ig*4]));
double gspikeg=1/(1+exp(-y[12+ig*4]/2));

f[12+ig*4]=(gnmdagg*nmdasig*(ENMDA-y[12+ig*4]) + gampag*ampasig*(EAMPA-y[12+ig*4])
	-glg1[ig]*(y[12+ig*4]+51)-gna*pow(minfgg,3)*y[14+ig*4]*(y[12+ig*4]-55)
    -gdrg*pow(y[13+ig*4],4)*(y[12+ig*4]+90) + dg*(gmf-y[12+ig*4]) + gachgaba*(EAch-y[12+ig*4]))*(1/Cgaba);
f[13+ig*4]=tng*(ang*(1-y[13+ig*4])-bng*y[13+ig*4]);
f[14+ig*4]=thg*(ahg*(1-y[14+ig*4])-bhg*y[14+ig*4]);
f[15+ig*4]=as*gspikeg*(1-y[15+ig*4])-bs*(1-gspikeg)*y[15+ig*4]; // GABA receptor activation;
}

}
void euler()
{
 int i=0;
 system();
  while (i<Neq){y[i]+=dt*f[i]; i++;}
} 

//nlhs - Number of expected output mxArrays
//plhs - Array of pointers to the expected output mxArrays
//nrhs - Number of input mxArrays
//prhs - Array of pointers to the input mxArrays. Do not modify any prhs values in your MEX-file. Changing the data in these read-only mxArrays can produce undesired side effects.

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	srand(1); //edit!   
//if(y) delete[] y;        
// y=new double[Neq];
mxArray *tempy = mxCreateDoubleMatrix(1,Neq,mxREAL);
y=mxGetPr(tempy);
    
//if(f) delete[] f;        
// f=new double[Neq];
mxArray *tempf = mxCreateDoubleMatrix(1,Neq,mxREAL);
f=mxGetPr(tempf);

double *input;
double *ggaba1, *gbarnmda1, *gampa1;

TT   = *(mxGetPr(prhs[0]));
input = mxGetPr(prhs[1]);
ggaba1 = (mxGetPr(prhs[2]));
gbarnmda1 = (mxGetPr(prhs[3]));
gampa1 = (mxGetPr(prhs[4]));

double *inputachda;
double *inputachgaba;
double *inputnic;
inputachda = mxGetPr(prhs[5]);
inputachgaba = mxGetPr(prhs[6]);
inputnic = mxGetPr(prhs[7]);
gbarachgaba = *(mxGetPr(prhs[8]));
gbarachda = *(mxGetPr(prhs[9]));
gL = *(mxGetPr(prhs[10]));
gbarh = *(mxGetPr(prhs[11]));

 int kk;
 tgl=0;
 N=(int)(TT/dt);
     
double *Vm;   
double *mas_Vgaba=new double[NG*N];
double *allgaba1;
double *Achcurrent;
double *nmdasig1;

plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[4] = mxCreateDoubleMatrix(1,NG*N,mxREAL);

    Vm = mxGetPr(plhs[0]);
   mas_Vgaba = mxGetPr(plhs[4]);
   allgaba1 = mxGetPr(plhs[1]);
   Achcurrent = mxGetPr(plhs[2]);
   nmdasig1 = mxGetPr(plhs[3]);



	 for(int i=0; i<NG;i++)
 { for(int j=0;j<N;j++) mas_Vgaba[i*N+j] = 0; }


    y[0]=-60.; y[1]=0.; y[2]=0; y[3]=50; y[4]=0.1; y[5]=0; y[6]=1.;y[7]=0; y[8]=0; y[9]=0; y[10]=0; y[11]=0.06;
 for (int ig=0; ig<NG; ig++){rngL[ig] = (double)rand()/RAND_MAX; y[12+ig*4]=-40+100*(rngL[ig]-0.5); glg1[ig]=(0.1+(0.1*(rngL[ig]-0.5)))*10; y[13+ig*4]=0; y[14+ig*4]=0; y[15+ig*4]=0.;}
	
   kk=0;
    for(int i=0; i<N; i++)
    {
        
        tgl+=dt;
      // inp=0;
       inp = input[i];
       inpachda=inputachda[i];
       inpachgaba=inputachgaba[i];
       inpnic=inputnic[i];
       ggaba=ggaba1[i];
       gbarnmda=gbarnmda1[i];
       gampa=gampa1[i];
        euler();
        Vm[kk]=y[0];
        allgaba1[i]=allgaba;
       Achcurrent[i]=gachda; //*(EAch-y[8]);
       nmdasig1[i]=y[1];
    
        
    for(int ig=0;ig<NG;ig++)
		{		
			mas_Vgaba[ig*N+kk]=y[12+ig*4];	
		}
         

        if(kk>N) break;
    kk++;
   }
   return;
} // end mexFunction()
