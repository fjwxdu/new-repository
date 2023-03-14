//*********************************************************************************************************************
//*********************************************************************************************************************
//*********************************************************************************************************************
//******																										*******
//******																										*******
//******		Light scattering of a spherical particle illuminated by	an elliptical Gaussian beam				*******
//******																										*******
//******														Email:jqshenk@163.com		2018/08/04			*******
//******																										*******
//*********************************************************************************************************************
//*********************************************************************************************************************
//******																										*******
//******	Negative time convention: exp(-iwt)																	*******
//******																										*******
//******	Normalized Associated Legendre function based on Ferrer's definition								*******
//******																										*******
//*********************************************************************************************************************
//*********************************************************************************************************************
//	references for BSC calculations:
//	1. Jianqi Shen, Xiaowei Jia, Haitao Yu, Compact formulation of the beam shape coefficients for elliptical Gaussian beam
//	   based on localized approximation, Journal of the Optical Society of America A 33(11), 2016: 2256-2263.
//	2. Wei Wang & Jianqi Shen: Beam shape coefficients calculation for an elliptical Gaussian beam with 1-dimensional quadrature
//	   and localized approximation methods, Journal of Quantitative Spectroscopy and Radiative Transfer 212, 2018: 139-148.
//	3. Jianqi Shen, Xiang Liu, Wei Wang, Haitao Yu: Calculation of light scattering of an elliptical Gaussian beam by a 
//	   spherical particle, Journal of the Optical Society of America A 35(8), 2018: 1288-1298.
//	references for Mie Coefficients:
//	1. Jianqi Shen & Xiaoshu Cai: Algorithm of Numerical Calculation on Lorentz Mie Theory, PIERS Online 1(6), 2005: 691-694.
//	Hints:
//	1. For large ratio of w0x and w0y, we recommend the ASD method for BSC calculation.
//	2. Any comments on this calculation or any problems/discussions are welcome from the authors, please contact by e-mail.

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
#include <time.h>
#include <dos.h>
#include <windows.h>
#include <complex.h>

const double pi=4.0*atan(1.0);
int NumberFuncCal=0,NumberHighDigit=0;

void Norm_FarField_Scattering_Intensity_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double CMIE[],double theta,double fai,double Isca[]);
void Norm_Incident_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double Einc[]);
void Norm_Internal_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double alpha,double rm,double im,double CMIE[],double Eint[]);
void Norm_Scatteri_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double alpha,double CMIE[],double Esca[]);
void Norm_Near_Field_Fast(int serm[],int sersp[],double VBSC[],int nminB,int nmaxB,int nmaxMie,double ruo,double theta,double fai,double alpha,double rm,double im, double CMIE[],double Enear[][3]);
void Sigma_m(int n,int serm[],int sersp[],double VBSC[],double NALFn0,double theta,double fai,double Sigmam[][4]);

void Modified_CMIE(int nmax,double rm,double im,double alpha,double CMIE[]);

int ReadDataBeam(int intway,double DataBeam[],int OrderBeam[]);
int Norm_BSC_Reading(int ctrl,char Pfile0[],char Pfile1[],int nmin,int nmax,int serm[],int sersp[],double VBSC[]);
int Norm_BSC_Writing(int method,char Pfile0[],char Pfile1[],int nmin,int nmax,double DataBeam[],int serm[]);
void Norm_BSC_LAM(int n,int m,double DataBeam[],double Anm[],double Bnm[]);
void Norm_BSC_Int(int intway,int n,int m,double DataBeam[],double Anm[],double Bnm[]);
void Integral_Q(int intway,int n,int m,double DataBeam[],double res[]);
int RombergTrapezoi_SubQ(int intway,int n,int m,double a,double b,double DataBeam[],double res[]);
int Integrand_Q1(int n,int m,double DataBeam[],double theta,double res[]);
void Integrand_Q0(int n,int m,double DataBeam[],double alpha,double res[]);
void FuncRelaModifiedBesselLAM(int m,double DblA[],double Beta[],double cosKsi[],double sinKsi[],double FImAnm[2],double FImBnm[2]);
int FuncRelaModifiedBessel(int m,double DblA[],double Beta[],double cosKsi[],double sinKsi[],double FIm[][2]);
void FuncXmBesselJB(int m,double DblA,double Beta,double Ksi,double Xmn1[],double Xmp1[]);
void FuncXmBesselJB1(int m,double DblA,double Beta,double Ksi,double Xmn1[],double Xmp1[]);

// BSCs of LG 
void Integraltm(int n,int m,double DataBeam[],double res[]);
int Romberg(int n,int m,double a,double b,double DataBeam[],double res[]);
void Integraltm0(double gamma,int n,int m,double DataBeam[],double res[]);
double def(int fx,int fy); 
int factorial(double x);
double Laguerre(double l_order,double p_order,double x);

void Norm_AngularFunc_Single(int n,int m,double theta,double NPAITAO[]);
void Norm_AngularFunc_COL(int n,double theta,double NPAI[],double NTAO[]);
void Norm_AssLegFunc_n0(int nmax,double theta,double NALFn0[]);
double Norm_AssLegFunc_Single(int n,int m, double theta);


void CMPLE_EqBesselI_3Order(int m,double rez,double imz,double resmn1[],double resm0[],double resmp1[]);
int CMPLE_EqBesselI(int nmax,double rez,double imz,double reBesselI[],double imBesselI[]);
int Real_EqBesselI(int nmax,double z,double BesselI[]);

void Real_BesselJ_ThreeOrder(int n, double z,double BslJ[]);
void Real_BesselJ_Downward(int nmax, double z,double BesselJ[]);
double Real_BesselJ_Downward_Single(int n, double z);

double Real_Riccati_BesselJ_Single(int n,double z);
void Real_Riccati_BesselJ(int nmax,double rez,double RB1[],double DRB1[]);
void Complex_EqRiccati_BesselJ(int nmax,double rez,double imz,double reRB1[],double imRB1[],double reDRB1[],double imDRB1[]);
void Real_Ricatti_BesselH(int nmax,double z,double reRB3[],double imRB3[],double reDRB3[],double imDRB3[]);
void LENTZ(int n,double rez,double imz,double res[]);

void CplxMul(double p,double q,double s,double t,double x[]);
void CplxDiv(double p,double q,double s,double t,double x[]);
void CplxPow(int p,double reR,double imR,double res[]);
void CplxLn(double a,double b,double res[]);
double Funcabs(double rez,double imz);




int main()
{
//	***************************************************************************************
	FILE *fp,*fp1,*fp2;
	int method,task[5];
	fp=fopen("Parameters for calculations.txt","r");
	fscanf(fp,"%d\n",&method);
	fscanf(fp,"%d\t%d\t%d\t%d\n",&task[0],&task[1],&task[2],&task[3]);
	fclose(fp);
//	***************************************************************************************
//	define the profile parameters:
	double DataBeam[20];
	int OrderBeam[3];
	ReadDataBeam(method,DataBeam,OrderBeam);
	int nminBeam=OrderBeam[0],nmaxBeam=OrderBeam[1];

//	***************************************************************************************
	char Pfile0[100],Pfile1[100],Pfile2[100],Pfile3[100],Pfile4[100],Pfile5[100],Pfile6[100];
	if(method==0)
	{
		strcpy(Pfile0,"ASD_BSC.txt");
		strcpy(Pfile1,"ASD_BSC_Record.txt");
		strcpy(Pfile2,"ASD_Inc_xy plane.csv");
		strcpy(Pfile3,"ASD_Far field scattering.csv");
		strcpy(Pfile4,"ASD_Inc_xz plane.csv");
		strcpy(Pfile5,"ASD_Sca_xz plane.csv");
		strcpy(Pfile6,"ASD_Tot_xz plane.csv");
	}
	else if(method==1)
	{
		strcpy(Pfile0,"Quad_BSC.txt");
		strcpy(Pfile1,"Quad_BSC_Record.txt");
		strcpy(Pfile2,"Quad_Inc_xy plane.csv");
		strcpy(Pfile3,"Quad_Far field scattering.csv");
		strcpy(Pfile4,"Quad_Inc_xz plane.csv");
		strcpy(Pfile5,"Quad_Sca_xz plane.csv");
		strcpy(Pfile6,"Quad_Tot_xz plane.csv");
	}
	else if(method==3)
	{
		strcpy(Pfile0,"piASD_BSC.txt");
		strcpy(Pfile1,"piASD_BSC_Record.txt");
		strcpy(Pfile2,"piASD_Inc_xy plane.csv");
		strcpy(Pfile3,"piASD_Far field scattering.csv");
		strcpy(Pfile4,"piASD_Inc_xz plane.csv");
		strcpy(Pfile5,"piASD_Sca_xz plane.csv");
		strcpy(Pfile6,"piASD_Tot_xz plane.csv");
	}
	else
	{
		strcpy(Pfile0,"LAM_BSC.txt");
		strcpy(Pfile1,"LAM_BSC_Record.txt");
		strcpy(Pfile2,"LAM_Inc_xy plane.csv");
		strcpy(Pfile3,"LAM_Far field scattering.csv");
		strcpy(Pfile4,"LAM_Inc_xz plane.csv");
		strcpy(Pfile5,"LAM_Sca_xz plane.csv");
		strcpy(Pfile6,"LAM_Tot_xz plane.csv");
	}

//	***************************************************************************************
//	***************************************************************************************
//	Calculate the BSCs:
//	***************************************************************************************
//	***************************************************************************************
	int *serm=(int*)malloc((nmaxBeam+5)*sizeof(int));
	int *sersp=(int*)malloc((nmaxBeam+5)*sizeof(int));
	double *VBSC=(double*)malloc(10*sizeof(double));
	int BSCreadok=0;
	if(task[0]==0)
	{
		fp=fopen(Pfile0,"r"); fp1=fopen(Pfile1,"r");
		if((!fp) || (!fp1))
		{
			task[0]=1;
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
			for(int warn=1;warn<=3;warn++) printf("\n\n\tBSCs are not ready, to be calculated soon...\n\n"),getche();
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
		}
		else
		{
			fclose(fp),fclose(fp1);
			BSCreadok=Norm_BSC_Reading(1,Pfile0,Pfile1,nminBeam,nmaxBeam,serm,sersp,VBSC);
			if(BSCreadok==0)
			{
				task[0]=1;
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
				for(int warn=1;warn<=3;warn++) printf("\n\n\tBSCs may be wrong, to be re-calculated soon...\n\n"),getche();
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
			}
		}
	}

	if(task[0]==1)
	{
		printf("\n\n******************************************************************\n\n");
		printf("\t����BSC��\tn = %d -> %d\n",nminBeam,nmaxBeam);
		printf("\n******************************************************************\n\n");
		clock_t start,finish;
		double Total_time;
		start=clock();
		Norm_BSC_Writing(method,Pfile0,Pfile1,nminBeam,nmaxBeam,DataBeam,serm);
		finish=clock();
		Total_time = (double)(finish-start) / CLOCKS_PER_SEC;

		fp=fopen(Pfile1,"a");
		if(method<=1 || method==3)
		{
			fprintf(fp,"%d\t%d\t%d\n",0,0,0);
			fprintf(fp,"\n\n\t CPU time for BSC calculation: %le(sec)\n",Total_time);
			printf("\n\n\tCPU time for BSC calculation: %le(sec)\n\n",Total_time);
			if(NumberHighDigit==0) fprintf(fp,"\tIntergand repeats %d times\n",NumberFuncCal);
			else fprintf(fp,"\tIntergand repeats %d*2147483647+%d times\n",NumberHighDigit,NumberFuncCal);
			printf("\tIntergand repeats %d times\n",NumberFuncCal);
		}
		else
		{
			fprintf(fp,"%d\t%d\t%d\n",0,0,0);
			fprintf(fp,"\n\n\tCPU time for BSC calculation: %le(sec)\n",Total_time);
			printf("\n\n\tCPU time for BSC calculation: %le(sec)\n\n",Total_time);
		}
		fclose(fp);
	}
	if(task[1]==0 && task[2]==0 && task[3]==0) return 0;

//	***************************************************************************************
//	***************************************************************************************
//	Read the BSCs:
//	***************************************************************************************
//	***************************************************************************************
	if(BSCreadok==0) BSCreadok=Norm_BSC_Reading(1,Pfile0,Pfile1,nminBeam,nmaxBeam,serm,sersp,VBSC);
	nmaxBeam=serm[0];
	VBSC=(double*)realloc(VBSC,(sersp[nmaxBeam+1]+5)*sizeof(double));
	Norm_BSC_Reading(2,Pfile0,Pfile1,nminBeam,nmaxBeam,serm,sersp,VBSC);


	double xf,yf,zf=0.0,xfmin,xfmax,yfmin,yfmax,zfmin,zfmax;
	double wavenumber0=DataBeam[1];
	double x0=DataBeam[14],y0=DataBeam[15],z0=DataBeam[16];

//	***************************************************************************************
//	***************************************************************************************
//	reconstruct indient beam in xy plane: zf=0
//	***************************************************************************************
//	***************************************************************************************
	if(task[1]==1)
	{
		printf("\n\n******************************************************************\n\n");
		printf("\t\treconstruct indient beam in xy plane");
		printf("\n******************************************************************\n\n");

		double scale=10.0;
		xfmax=1.0+int(fabs(x0)+scale*DataBeam[17]);
		yfmax=1.0+int(fabs(y0)+scale*DataBeam[18]);
		if(xfmax<yfmax) xfmax=yfmax;
		else yfmax=xfmax;
		xfmin=-xfmax,yfmin=-yfmax;

		double ruo,theta=0.0,fai=0.0,Einc[5],Vfield;
		int ix,iy,px=200,py=200;
		double stepx=(xfmax-xfmin)/px,stepy=(yfmax-yfmin)/py;

		fp=fopen(Pfile2,"w");
		if(method==0) fprintf(fp,"ASD_inc");
		else if(method==1) fprintf(fp,"Int_inc");
		else if(method==3) fprintf(fp,"piASD_inc");
		else  fprintf(fp,"LAM_Inc");
		for(ix=0;ix<=px;ix++) fprintf(fp,"\t%lf",xfmin+ix*stepx);
		fprintf(fp,"\n");
		for(iy=0;iy<=py;iy++)
		{
			yf=yfmin+iy*stepy;
			fprintf(fp,"%lf",yf);
			for(Vfield=0.0,ix=0;ix<=px;ix++)
			{
				xf=xfmin+ix*stepx;
				double thetaold=theta,faiold=fai;
				double radial=sqrt(xf*xf+yf*yf+zf*zf);
				ruo=radial*wavenumber0;
				theta=(radial>0.0)? acos(zf/radial):thetaold;
				double radialxy=sqrt(xf*xf+yf*yf);
				fai=(radialxy>0.0)? ((yf>=0.0)? acos(xf/radialxy):-acos(xf/radialxy)):faiold;
				if(nminBeam>1 && ruo==0.0) ruo=1e-5;

				Norm_Incident_Field_Fast(serm,sersp,VBSC,nminBeam,nmaxBeam,ruo,theta,fai,Einc);
				printf("\r\tiy= %3d / %3d\tix= %3d / %3d",iy,py,ix,px);
				
				fprintf(fp,"\t%le",Einc[0]);
				if(Vfield<Einc[0]) Vfield=Einc[0];
			}
			fprintf(fp,"\n");
			fclose(fp);fp=fopen(Pfile2,"a");
			printf("\r\tiy= %3d / %d\tyf= %.3f\tMax = %.1e\n",iy,py,yf,Vfield);
		}
		fclose(fp);
	}


//	***************************************************************************************
//	***************************************************************************************
//	Calculate the scattering
//	***************************************************************************************
//	***************************************************************************************
	if(task[2]!=0 || task[3]!=0)
	{
//	Calculate the Mie Coefficients
		int theoryctrl=0,Debyepmin=-1,Debyepmax=-1;
		double rm,im,diaparticle;
		fp=fopen("Parameters for particle.txt","r");
		fscanf(fp,"%lf\n",&diaparticle);
		fscanf(fp,"%lf\t%lf\n",&rm,&im);
		fscanf(fp,"%d\n",&theoryctrl);
		if(theoryctrl>0) fscanf(fp,"%d\t%d\n",&Debyepmin,&Debyepmax);
		fclose(fp);
		double alpha=0.5*diaparticle*wavenumber0;
		int nmaxMie=int(alpha+7.5*pow(alpha,0.4))+20;
		int nmax1=nmaxMie*8+5;
		double *CMIE=(double*)calloc(nmax1,sizeof(double));
		Modified_CMIE(nmaxMie,rm,im,alpha,CMIE);
		int nmaxSca=(nmaxBeam<nmaxMie)? nmaxBeam:nmaxMie;

//	Calculate the far field scattering in xz plane
		if(task[2]!=0)
		{
			printf("\n\n******************************************************************\n\n");
			printf("\t\tCalculate the far field scattering");
			printf("\n******************************************************************\n\n");
			double Isca[3],fai,theta,tamax=pi;
			fp=fopen(Pfile3,"w");
			if(method==0) fprintf(fp,"angle\tI1_ASD\tI2_ASD\n");
			else if(method==1) fprintf(fp,"angle\tI1_Int\tI2_Int\n");
			else fprintf(fp,"angle\tI1_LAM\tI2_LAM\n");
			for(int ita=-1800;ita<=1800;ita++)
			{
				fai=(ita<0)? pi:0;
				theta=tamax*fabs(ita)/1800.0;
				double thetap=ita*0.1;
				Norm_FarField_Scattering_Intensity_Fast(serm,sersp,VBSC,nminBeam,nmaxSca,CMIE,theta,fai,Isca);
				printf("\t%lf\t%le\t%le\n",thetap,Isca[1],Isca[2]);
				fprintf(fp,"%lf\t%le\t%le\n",thetap,Isca[1],Isca[2]);
			}
			fclose(fp);
		}

//	Calculate the near field in xz plane
		if(task[3]!=0)
		{
			double scale=1.5;
			printf("\n\n******************************************************************\n\n");
			printf("\t\tCalculate the near field in xz plane");
			printf("\n******************************************************************\n\n");
			xfmax=zfmax=1.0+int(scale*diaparticle),xfmin=zfmin=-xfmax,yf=0.0;
			double ruo,theta=0,fai=0,Enear[5][3],Vfield;
			int ix,iz,px=200,pz=200;
			double stepx=(xfmax-xfmin)/px,stepz=(zfmax-zfmin)/pz;
			fp=fopen(Pfile4,"w"); fp1=fopen(Pfile5,"w"); fp2=fopen(Pfile6,"w");
			if(method==0) fprintf(fp,"ADS_Inc"),fprintf(fp1,"ADS_Sca"),fprintf(fp2,"ADS_Tot");
			else if(method==1) fprintf(fp,"Int_Inc"),fprintf(fp1,"Int_Sca"),fprintf(fp2,"Int_Tot");
			else fprintf(fp,"LAM_Inc"),fprintf(fp1,"LAM_Sca"),fprintf(fp2,"LAM_Tot");
			for(ix=0;ix<=px;ix++) fprintf(fp,"\t%lf",xfmin+ix*stepx),fprintf(fp1,"\t%lf",xfmin+ix*stepx),fprintf(fp2,"\t%lf",xfmin+ix*stepx);
			fprintf(fp,"\n"),fprintf(fp1,"\n"),fprintf(fp2,"\n");
			for(iz=0;iz<=pz;iz++)
			{
				zf=zfmin+iz*stepz;
				fprintf(fp,"%lf",zf); fprintf(fp1,"%lf",zf); fprintf(fp2,"%lf",zf);
				for(Vfield=0.0,ix=0;ix<=px;ix++)
				{
					xf=xfmin+ix*stepx;
					double thetaold=theta,faiold=fai;
					double radial=sqrt(xf*xf+yf*yf+zf*zf);
					ruo=radial*wavenumber0;
					theta=(radial>0.0)? acos(zf/radial):thetaold;
					double radialxy=sqrt(xf*xf+yf*yf);
					fai=(radialxy>0.0)? ((yf>=0.0)? acos(xf/radialxy):-acos(xf/radialxy)):faiold;
					if(nminBeam>1 && ruo==0.0) ruo=1e-5;
					Norm_Near_Field_Fast(serm,sersp,VBSC,nminBeam,nmaxBeam,nmaxMie,ruo,theta,fai,alpha,rm,im,CMIE,Enear);
					printf("\r\tiz=%3d / %d\tix=%3d / %d",iz,pz,ix,px);
					fprintf(fp,"\t%le",Enear[0][0]); fprintf(fp1,"\t%le",Enear[0][1]); fprintf(fp2,"\t%le",Enear[0][2]);
					if(Vfield<Enear[0][0]) Vfield=Enear[0][0];
					if(Vfield<Enear[0][1]) Vfield=Enear[0][1];
					if(Vfield<Enear[0][2]) Vfield=Enear[0][2];
				}
				fprintf(fp,"\n"),fprintf(fp1,"\n"),fprintf(fp2,"\n");
				fclose(fp),fclose(fp1),fclose(fp2);
				fp=fopen(Pfile4,"a"),fp1=fopen(Pfile5,"a"),fp2=fopen(Pfile6,"a");
				printf("\r\tiz= %3d / %3d\tzf= %.3f\tMax = %.3e\n",iz,pz,zf,Vfield);
			}
			fclose(fp),fclose(fp1),fclose(fp2);
		}
		free(CMIE); CMIE=NULL;
	}
	free(serm),free(sersp),free(VBSC); VBSC=NULL,serm=NULL,sersp=NULL;
}



//	***************************************************************************************
//	***************************************************************************************
//	Parameters for beam profile
//	***************************************************************************************
//	***************************************************************************************
int ReadDataBeam(int intway,double DataBeam[],int OrderBeam[])
{
	double w0x,w0y,wavel0,x0,y0,z0;
	double l_order,p_order;
	FILE *fp;
	fp=fopen("Parameters of incident beam.txt","r");
	fscanf(fp,"%lf\t%lf\t\t%lf\n",&w0x,&w0y,&wavel0);
	fscanf(fp,"%lf\t%lf\t%lf\n",&x0,&y0,&z0);
	if(intway==3){
		fscanf(fp,"%lf\t%lf\n",&l_order,&p_order);
	}
	fclose(fp);
//	***************************************************************************************
	double wavenumber0=2.0*pi/wavel0;
	double X0=wavenumber0*x0,Y0=wavenumber0*y0,Z0=wavenumber0*z0;
	double sx=1.0/wavenumber0/w0x,sy=1.0/wavenumber0/w0y;
//	***************************************************************************************
//	nminBeam��nmaxBeam��
	double zRx=pi*w0x*w0x/wavel0,zRy=pi*w0y*w0y/wavel0;
	double wzx=w0x*sqrt(1.0+pow(z0/zRx,2));
	double wzy=w0y*sqrt(1.0+pow(z0/zRy,2));
	double Wzx=wzx*wavenumber0,Wzy=wzy*wavenumber0;
	double R0=sqrt(X0*X0+Y0*Y0);
	double wratio=sqrt(10.0*log(10.0));
	double wwx=Wzx*wratio,wwy=Wzy*wratio,Rmax=0.0;
	for(int ialpha=0;ialpha<720;ialpha++)
	{
		double alpha=ialpha*pi/360.0;
		double Ralpha=sqrt(pow(X0+wwx*cos(alpha),2)+pow(Y0+wwy*sin(alpha),2));
		if(Rmax<Ralpha) Rmax=Ralpha;
	}
	int nmaxBeam=int(Rmax)+1;
	int nminBeam=1;
	int nmidBeam=int(R0+0.5);

//	***************************************************************************************
//	Keywords: Beam-Profile
//	DataBeam[0]		ruon=n+0.5
//	DataBeam[1]		wavenumber0
//	DataBeam[2->?]
	DataBeam[1]=wavenumber0,DataBeam[2]=sx,DataBeam[3]=sy,DataBeam[4]=X0,DataBeam[5]=Y0,DataBeam[6]=Z0,DataBeam[7]=Wzx,DataBeam[8]=Wzy;
	if(intway==3){
		DataBeam[9]=l_order;DataBeam[10]=p_order;
	}
	DataBeam[11]=wavel0,DataBeam[12]=w0x,DataBeam[13]=w0y,DataBeam[14]=x0,DataBeam[15]=y0,DataBeam[16]=z0,DataBeam[17]=wzx,DataBeam[18]=wzy;
	OrderBeam[0]=nminBeam,OrderBeam[1]=nmaxBeam,OrderBeam[2]=nmidBeam;
	return 1;
}



//	***************************************************************************************
//	***************************************************************************************
//	Far field scattering
//	***************************************************************************************
//	***************************************************************************************
void Norm_FarField_Scattering_Intensity_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double CMIE[],double theta,double fai,double Isca[])
{
	if(nmin>nmax) {Isca[1]=Isca[2]=0.0;return;}
//	***************************************************************************************
	double *NPAI=(double*)malloc((nmax+5)*sizeof(double)),*NTAO=(double*)malloc((nmax+5)*sizeof(double));

	double tmp1[2],tmp2[2],reS1=0.0,imS1=0.0,reS2=0.0,imS2=0.0;
	for(int n=nmin;n<=nmax;n++)
	{
		Norm_AngularFunc_COL(n,theta,NPAI,NTAO);
		double *BSCn0=&VBSC[sersp[n]];
		double reS11=0.0,imS11=0.0,reS12=NTAO[0]*BSCn0[3],imS12=NTAO[0]*BSCn0[4];
		double reS21=NTAO[0]*BSCn0[1],imS21=NTAO[0]*BSCn0[2],reS22=0.0,imS22=0.0;
		for(int m=1;m<=serm[n];m++)
		{
			double *BSCnm=&VBSC[sersp[n]+m*8];
			double cosmfai=cos(m*fai),sinmfai=sin(m*fai);
			reS11+=NPAI[m]*((BSCnm[2]-BSCnm[6])*cosmfai+(BSCnm[1]+BSCnm[5])*sinmfai);
			imS11+=NPAI[m]*((BSCnm[5]-BSCnm[1])*cosmfai+(BSCnm[2]+BSCnm[6])*sinmfai);
			reS12+=NTAO[m]*((BSCnm[3]+BSCnm[7])*cosmfai+(BSCnm[8]-BSCnm[4])*sinmfai);
			imS12+=NTAO[m]*((BSCnm[4]+BSCnm[8])*cosmfai+(BSCnm[3]-BSCnm[7])*sinmfai);

			reS21+=NTAO[m]*((BSCnm[1]+BSCnm[5])*cosmfai+(BSCnm[6]-BSCnm[2])*sinmfai);
			imS21+=NTAO[m]*((BSCnm[2]+BSCnm[6])*cosmfai+(BSCnm[1]-BSCnm[5])*sinmfai);
			reS22+=NPAI[m]*((BSCnm[8]-BSCnm[4])*cosmfai-(BSCnm[3]+BSCnm[7])*sinmfai);
			imS22+=NPAI[m]*((BSCnm[3]-BSCnm[7])*cosmfai-(BSCnm[4]+BSCnm[8])*sinmfai);
		}
//		***********************************************************************************
		double *CMIEn=&CMIE[8*(n-1)];
		double paran=(n+0.5)/n/(n+1.0);
		CplxMul(CMIEn[1],CMIEn[2],reS11,imS11,tmp1);
		CplxMul(CMIEn[3],CMIEn[4],reS12,imS12,tmp2);
		reS1+=(tmp1[0]+tmp2[0])*paran,imS1+=(tmp1[1]+tmp2[1])*paran;

		CplxMul(CMIEn[1],CMIEn[2],reS21,imS21,tmp1);
		CplxMul(CMIEn[3],CMIEn[4],reS22,imS22,tmp2);
		reS2+=(tmp1[0]+tmp2[0])*paran,imS2+=(tmp1[1]+tmp2[1])*paran;
//		***********************************************************************************
	}
	free(NPAI),free(NTAO); NPAI=NTAO=NULL;
	Isca[1]=reS1*reS1+imS1*imS1,Isca[2]=reS2*reS2+imS2*imS2;
	return;
}




//	***************************************************************************************
//	***************************************************************************************
//	Field inside the spherical particle
//	***************************************************************************************
//	***************************************************************************************
void Norm_Internal_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double alpha,double rm,double im,double CMIE[],double Eint[])
{
	if(nmin>nmax || ruo>alpha) {Eint[0]=Eint[1]=Eint[2]=Eint[3]=0.0;return;}
//	***************************************************************************************
	if(ruo>0.0)
	{
		int nmax1=nmax+5;
		double *NALFn0=(double*)malloc(nmax1*sizeof(double));
		double *reRB1=(double*)malloc(nmax1*sizeof(double)),*imRB1=(double*)malloc(nmax1*sizeof(double));
		double *reDRB1=(double*)malloc(nmax1*sizeof(double)),*imDRB1=(double*)malloc(nmax1*sizeof(double));

		double reruosp=rm*ruo,imruosp=im*ruo;
		Complex_EqRiccati_BesselJ(nmax,reruosp,imruosp,reRB1,imRB1,reDRB1,imDRB1);
		Norm_AssLegFunc_n0(nmax,theta,NALFn0);

		double Sigmam[3][4],reSR=0.0,imSR=0.0,reST=0.0,imST=0.0,reSF=0.0,imSF=0.0;
		for(int n=nmin;n<=nmax;n++)
		{
			Sigma_m(n,serm,sersp,VBSC,NALFn0[n],theta,fai,Sigmam);
			double repwr,impwr,repwtf,impwtf;
			switch(n%4)
			{
				case 0: repwr=n+0.50,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0; break;
				case 1: impwr=n+0.50,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0; break;
				case 2: repwr=-n-0.5,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0; break;
				case 3: impwr=-n-0.5,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0; break;
			}
//		***********************************************************************************
			double tmp[2],tmp0[2],tmp1[2],tmp2[2],tmp3[2],tmp4[2];
			double *CMIEn=&CMIE[8*(n-1)];
			CplxMul(CMIEn[5],CMIEn[6], reRB1[n], imRB1[n],tmp0);	//	Cn* RB1
			CplxMul(CMIEn[5],CMIEn[6],reDRB1[n],imDRB1[n],tmp1);	//	Cn*DRB1
			CplxMul(CMIEn[7],CMIEn[8], reRB1[n], imRB1[n],tmp2);	//	Dn* RB1
//		***********************************************************************************
			CplxMul(repwr,impwr,tmp0[0],tmp0[1],tmp);
			CplxMul(tmp[0],tmp[1],Sigmam[0][0],Sigmam[0][1],tmp);
			reSR-=tmp[1],imSR+=tmp[0];

			CplxMul(tmp1[0],tmp1[1],Sigmam[1][0],Sigmam[1][1],tmp3);
			CplxMul(tmp2[0],tmp2[1],Sigmam[1][2],Sigmam[1][3],tmp4);
			CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
			reST-=tmp[1],imST+=tmp[0];

			CplxMul(tmp2[0],tmp2[1],Sigmam[2][2],Sigmam[2][3],tmp3);
			CplxMul(tmp1[0],tmp1[1],Sigmam[2][0],Sigmam[2][1],tmp4);
			CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
			reSF-=tmp[0],imSF-=tmp[1];
		}
		double ruosquare=reruosp*reruosp+imruosp*imruosp;
		Eint[1]=(reSR*reSR+imSR*imSR)/ruosquare/ruosquare;		//	radial component
		Eint[2]=(reST*reST+imST*imST)/ruosquare;				//	theta
		Eint[3]=(reSF*reSF+imSF*imSF)/ruosquare;				//	fai
		Eint[0]=Eint[1]+Eint[2]+Eint[3];						//	total field in square

		free(NALFn0),free(reRB1),free(imRB1),free(reDRB1),free(imDRB1); NALFn0=reRB1=imRB1=reDRB1=imDRB1=NULL;
	}
//	***************************************************************************************
	else
	{
		double Mie_c1square=CMIE[5]*CMIE[5]+CMIE[6]*CMIE[6];
		Eint[0]=Mie_c1square*(VBSC[13]*VBSC[13]+VBSC[14]*VBSC[14]+VBSC[9]*VBSC[9]+VBSC[10]*VBSC[10]+VBSC[1]*VBSC[1]+VBSC[2]*VBSC[2])*0.375;
	}
	double Modifyfactor=exp(2.0*im*(ruo-alpha));
	for(int i=0;i<=3;i++) Eint[i]=sqrt(Eint[i]*Modifyfactor);
	return;
}





//	***************************************************************************************
//	***************************************************************************************
//	scattered field
//	***************************************************************************************
//	***************************************************************************************
void Norm_Scatteri_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double alpha,double CMIE[],double Esca[])
{
	if(nmin>nmax || ruo<=alpha) {Esca[0]=Esca[1]=Esca[2]=Esca[3]=0.0;return;}
//	***************************************************************************************
	int nmax1=nmax+5;
	double *NALFn0=(double*)malloc(nmax1*sizeof(double));
	double *reRB3=(double*)malloc(nmax1*sizeof(double)),*imRB3=(double*)malloc(nmax1*sizeof(double));
	double *reDRB3=(double*)malloc(nmax1*sizeof(double)),*imDRB3=(double*)malloc(nmax1*sizeof(double));

	Norm_AssLegFunc_n0(nmax,theta,NALFn0);
	Real_Ricatti_BesselH(nmax,ruo,reRB3,imRB3,reDRB3,imDRB3);

//	***************************************************************************************
	double Sigmam[3][4],tmp[2],reSR=0.0,imSR=0.0,reST=0.0,imST=0.0,reSF=0.0,imSF=0.0;
	for(int n=nmin;n<=nmax;n++)
	{
		Sigma_m(n,serm,sersp,VBSC,NALFn0[n],theta,fai,Sigmam);
		double repwr,impwr,repwtf,impwtf;
		switch(n%4)
		{
			case 0: repwr=n+0.50,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
			case 1: impwr=n+0.50,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
			case 2: repwr=-n-0.5,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
			case 3: impwr=-n-0.5,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
		}
//		***********************************************************************************
		double tmp0[2],tmp1[2],tmp2[2],tmp3[2],tmp4[2];
		double *CMIEn=&CMIE[8*(n-1)];
		CplxMul(CMIEn[1],CMIEn[2], reRB3[n], imRB3[n],tmp0);	//	an* RB3
		CplxMul(CMIEn[1],CMIEn[2],reDRB3[n],imDRB3[n],tmp1);	//	an*DRB3
		CplxMul(CMIEn[3],CMIEn[4], reRB3[n], imRB3[n],tmp2);	//	bn* RB3
//		***********************************************************************************

		CplxMul(repwr,impwr,tmp0[0],tmp0[1],tmp);
		CplxMul(tmp[0],tmp[1],Sigmam[0][0],Sigmam[0][1],tmp);
		reSR+=tmp[1],imSR-=tmp[0];

		CplxMul(tmp1[0],tmp1[1],Sigmam[1][0],Sigmam[1][1],tmp3);
		CplxMul(tmp2[0],tmp2[1],Sigmam[1][2],Sigmam[1][3],tmp4);
		CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
		reST+=tmp[1],imST-=tmp[0];

		CplxMul(tmp2[0],tmp2[1],Sigmam[2][2],Sigmam[2][3],tmp3);
		CplxMul(tmp1[0],tmp1[1],Sigmam[2][0],Sigmam[2][1],tmp4);
		CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
		reSF+=tmp[0],imSF+=tmp[1];
	}
	double ruosquare=ruo*ruo;
	Esca[1]=(reSR*reSR+imSR*imSR)/ruosquare/ruosquare;
	Esca[2]=(reST*reST+imST*imST)/ruosquare;
	Esca[3]=(reSF*reSF+imSF*imSF)/ruosquare;
	Esca[0]=Esca[1]+Esca[2]+Esca[3];
	free(NALFn0),free(reRB3),free(imRB3),free(reDRB3),free(imDRB3); NALFn0=reRB3=imRB3=reDRB3=imDRB3=NULL;
	for(int i=0;i<=3;i++) Esca[i]=sqrt(Esca[i]);
	return;
}




//	***************************************************************************************
//	***************************************************************************************
//	field calculation for incident beam
//	***************************************************************************************
//	***************************************************************************************
void Norm_Incident_Field_Fast(int serm[],int sersp[],double VBSC[],int nmin,int nmax,double ruo,double theta,double fai,double Einc[])
{
	if(nmin>nmax) {Einc[0]=Einc[1]=Einc[2]=Einc[3]=0.0;return;}
//	***************************************************************************************
	if(ruo>0.0)
	{
		int nmax1=nmax+5;
		double *NALFn0=(double*)malloc(nmax1*sizeof(double));
		double *RB1=(double*)malloc(nmax1*sizeof(double)),*DRB1=(double*)malloc(nmax1*sizeof(double));
		Norm_AssLegFunc_n0(nmax,theta,NALFn0);
		Real_Riccati_BesselJ(nmax,ruo,RB1,DRB1);

		double Sigmam[3][4],tmp[2],reSR=0.0,imSR=0.0,reST=0.0,imST=0.0,reSF=0.0,imSF=0.0;
		for(int n=nmin;n<=nmax;n++)
		{
			Sigma_m(n,serm,sersp,VBSC,NALFn0[n],theta,fai,Sigmam);
			double repwr,impwr,repwtf,impwtf;
			switch(n%4)
			{
				case 0: repwr=n+0.50,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
				case 1: impwr=n+0.50,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
				case 2: repwr=-n-0.5,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
				case 3: impwr=-n-0.5,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
			}
//		***********************************************************************************
			CplxMul(repwr,impwr,Sigmam[0][0]*RB1[n],Sigmam[0][1]*RB1[n],tmp);
			reSR-=tmp[1],imSR+=tmp[0];

			CplxMul(repwtf,impwtf,Sigmam[1][0]*DRB1[n]-Sigmam[1][2]*RB1[n],Sigmam[1][1]*DRB1[n]-Sigmam[1][3]*RB1[n],tmp);
			reST-=tmp[1],imST+=tmp[0];

			CplxMul(repwtf,impwtf,Sigmam[2][2]*RB1[n]-Sigmam[2][0]*DRB1[n],Sigmam[2][3]*RB1[n]-Sigmam[2][1]*DRB1[n],tmp);
			reSF-=tmp[0],imSF-=tmp[1];
		}
		double ruosquare=ruo*ruo;
		Einc[1]=(reSR*reSR+imSR*imSR)/ruosquare/ruosquare;
		Einc[2]=(reST*reST+imST*imST)/ruosquare;
		Einc[3]=(reSF*reSF+imSF*imSF)/ruosquare;
		Einc[0]=Einc[1]+Einc[2]+Einc[3];
//		Einc[0]=sqrt(Einc[0]); 
//		Einc[0]=sqrt(Einc[1])*sin(theta)*cos(fai)+sqrt(Einc[2])*cos(theta)*cos(fai)-sqrt(Einc[3])*sin(fai);
		
		free(NALFn0),free(RB1),free(DRB1); NALFn0=RB1=DRB1=NULL;
	}
//	***************************************************************************************
	else Einc[0]=(VBSC[13]*VBSC[13]+VBSC[14]*VBSC[14]+VBSC[9]*VBSC[9]+VBSC[10]*VBSC[10]+VBSC[1]*VBSC[1]+VBSC[2]*VBSC[2])*0.375;
	for(int i=0;i<=3;i++) Einc[i]=sqrt(Einc[i]);
	return;
}



//	***************************************************************************************
//	***************************************************************************************
//	to calculate the incident field, scattered/internal field and the total field 
//	Inc: Enear[0][0]
//	Sca: Enear[0][1]			ruo<=alpha, internal		ruo>alpha, scattered 
//	Tot: Enear[0][2]			ruo<=alpha, internal		ruo>alpha, Tot=Inc+Sca
//	***************************************************************************************
//	***************************************************************************************
void Norm_Near_Field_Fast(int serm[],int sersp[],double VBSC[],int nminB,int nmaxB,int nmaxMie,double ruo,double theta,double fai,double alpha,double rm,double im, double CMIE[],double Enear[][3])
{
	int i,j;
	if(nminB>nmaxB) {for(i=0;i<4;i++) {for(j=0;j<3;j++) Enear[i][j]=0.0;} return;}
//	***************************************************************************************
//		��������ԭ������ĳ�
//	***************************************************************************************
	if(ruo>0.0)
	{
		int nmax1=nmaxB+5;
		double *NALFn0=(double*)calloc(nmax1,sizeof(double));
		double *RBInc=(double*)calloc(nmax1,sizeof(double)),*DRBInc=(double*)calloc(nmax1,sizeof(double));
		double *reRBSca=(double*)calloc(nmax1,sizeof(double)),*imRBSca=(double*)calloc(nmax1,sizeof(double));
		double *reDRBSca=(double*)calloc(nmax1,sizeof(double)),*imDRBSca=(double*)calloc(nmax1,sizeof(double));

		Real_Riccati_BesselJ(nmaxB,ruo,RBInc,DRBInc);
		if(ruo<=alpha) Complex_EqRiccati_BesselJ(nmaxB,rm*ruo,im*ruo,reRBSca,imRBSca,reDRBSca,imDRBSca);
		else Real_Ricatti_BesselH(nmaxB,ruo,reRBSca,imRBSca,reDRBSca,imDRBSca);
		Norm_AssLegFunc_n0(nmaxB,theta,NALFn0);

		double Sigmam[3][4],tmp[2],tmp0[2],tmp1[2],tmp2[2],tmp3[2],tmp4[2];
		double reSR[3]={0,0,0},imSR[3]={0,0,0},reST[3]={0,0,0},imST[3]={0,0,0},reSF[3]={0,0,0},imSF[3]={0,0,0};
		for(int n=nminB;n<=nmaxB;n++)
		{
			Sigma_m(n,serm,sersp,VBSC,NALFn0[n],theta,fai,Sigmam);
			double repwr,impwr,repwtf,impwtf;
			switch(n%4)
			{
				case 0: repwr=n+0.50,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
				case 1: impwr=n+0.50,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
				case 2: repwr=-n-0.5,repwtf=repwr/n/(n+1.0),impwr=impwtf=0.0;break;
				case 3: impwr=-n-0.5,impwtf=impwr/n/(n+1.0),repwr=repwtf=0.0;break;
			}

//		***********************************************************************************
			CplxMul(repwr,impwr,Sigmam[0][0]*RBInc[n],Sigmam[0][1]*RBInc[n],tmp);
			reSR[0]-=tmp[1],imSR[0]+=tmp[0];

			CplxMul(repwtf,impwtf,Sigmam[1][0]*DRBInc[n]-Sigmam[1][2]*RBInc[n],Sigmam[1][1]*DRBInc[n]-Sigmam[1][3]*RBInc[n],tmp);
			reST[0]-=tmp[1],imST[0]+=tmp[0];

			CplxMul(repwtf,impwtf,Sigmam[2][2]*RBInc[n]-Sigmam[2][0]*DRBInc[n],Sigmam[2][3]*RBInc[n]-Sigmam[2][1]*DRBInc[n],tmp);
			reSF[0]-=tmp[0],imSF[0]-=tmp[1];
//		***********************************************************************************
			if(n<=nmaxMie)
			{
				double *CMIEn=&CMIE[8*(n-1)];
				if(ruo<=alpha)
				{
					CplxMul(CMIEn[5],CMIEn[6], reRBSca[n], imRBSca[n],tmp0);	//	Cn* RBN
					CplxMul(CMIEn[5],CMIEn[6],reDRBSca[n],imDRBSca[n],tmp1);	//	Cn*DRBN
					CplxMul(CMIEn[7],CMIEn[8], reRBSca[n], imRBSca[n],tmp2);	//	Dn* RBN
					CplxMul(repwr,impwr,tmp0[0],tmp0[1],tmp);
					CplxMul(tmp[0],tmp[1],Sigmam[0][0],Sigmam[0][1],tmp);
					reSR[1]-=tmp[1],imSR[1]+=tmp[0];

					CplxMul(tmp1[0],tmp1[1],Sigmam[1][0],Sigmam[1][1],tmp3);
					CplxMul(tmp2[0],tmp2[1],Sigmam[1][2],Sigmam[1][3],tmp4);
					CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
					reST[1]-=tmp[1],imST[1]+=tmp[0];

					CplxMul(tmp2[0],tmp2[1],Sigmam[2][2],Sigmam[2][3],tmp3);
					CplxMul(tmp1[0],tmp1[1],Sigmam[2][0],Sigmam[2][1],tmp4);
					CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
					reSF[1]-=tmp[0],imSF[1]-=tmp[1];
				}
				else
				{
					CplxMul(CMIEn[1],CMIEn[2], reRBSca[n], imRBSca[n],tmp0);	//	an* RB3
					CplxMul(CMIEn[1],CMIEn[2],reDRBSca[n],imDRBSca[n],tmp1);	//	an*DRB3
					CplxMul(CMIEn[3],CMIEn[4], reRBSca[n], imRBSca[n],tmp2);	//	bn* RB3

					CplxMul(repwr,impwr,tmp0[0],tmp0[1],tmp);
					CplxMul(tmp[0],tmp[1],Sigmam[0][0],Sigmam[0][1],tmp);
					reSR[1]+=tmp[1],imSR[1]-=tmp[0];

					CplxMul(tmp1[0],tmp1[1],Sigmam[1][0],Sigmam[1][1],tmp3);
					CplxMul(tmp2[0],tmp2[1],Sigmam[1][2],Sigmam[1][3],tmp4);
					CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
					reST[1]+=tmp[1],imST[1]-=tmp[0];

					CplxMul(tmp2[0],tmp2[1],Sigmam[2][2],Sigmam[2][3],tmp3);
					CplxMul(tmp1[0],tmp1[1],Sigmam[2][0],Sigmam[2][1],tmp4);
					CplxMul(repwtf,impwtf,tmp3[0]-tmp4[0],tmp3[1]-tmp4[1],tmp);
					reSF[1]+=tmp[0],imSF[1]+=tmp[1];
				}
			}
		}
//		***************************************************************************************
		if(ruo>alpha)
		{
			reSR[2]=reSR[1]+reSR[0],reST[2]=reST[1]+reST[0],reSF[2]=reSF[1]+reSF[0];
			imSR[2]=imSR[1]+imSR[0],imST[2]=imST[1]+imST[0],imSF[2]=imSF[1]+imSF[0];
		}
		else reSR[2]=reSR[1],reST[2]=reST[1],reSF[2]=reSF[1],imSR[2]=imSR[1],imST[2]=imST[1],imSF[2]=imSF[1];
//		***************************************************************************************
		for(j=0;j<=2;j++)
		{
			double ruosquare=ruo*ruo;
			if(j!=0 && ruo<=alpha) ruosquare*=(rm*rm+im*im);
			Enear[1][j]=(reSR[j]*reSR[j]+imSR[j]*imSR[j])/ruosquare/ruosquare;
			Enear[2][j]=(reST[j]*reST[j]+imST[j]*imST[j])/ruosquare;
			Enear[3][j]=(reSF[j]*reSF[j]+imSF[j]*imSF[j])/ruosquare;
			Enear[0][j]=Enear[1][j]+Enear[2][j]+Enear[3][j];
		}
		free(NALFn0),free(RBInc),free(DRBInc); NALFn0=RBInc=DRBInc=NULL;
		free(reRBSca),free(imRBSca),free(reDRBSca),free(imDRBSca);reRBSca=imRBSca=reDRBSca=imDRBSca=NULL;
	}
//	***************************************************************************************
	else
	{
		Enear[0][0]=(VBSC[1]*VBSC[1]+VBSC[2]*VBSC[2]+VBSC[9]*VBSC[9]+VBSC[10]*VBSC[10]+VBSC[13]*VBSC[13]+VBSC[14]*VBSC[14])*0.375;
		double Mie_c1square=CMIE[5]*CMIE[5]+CMIE[6]*CMIE[6];
		Enear[0][2]=Enear[0][1]=Enear[0][0]*Mie_c1square;
	}

	if(ruo<=alpha)
	{
		double Modifyfactor=exp(2.0*im*(ruo-alpha));
		for(i=0;i<=3;i++) Enear[i][1]*=Modifyfactor,Enear[i][2]=Enear[i][1];
	}
	for(i=0;i<=3;i++) {for(j=0;j<=2;j++) Enear[i][j]=sqrt(Enear[i][j]);}
	return;
}




//	**************************************************************************************************
//	**************************************************************************************************
//	Sum for azimuth mode m
//	**************************************************************************************************
//	**************************************************************************************************
void Sigma_m(int n,int serm[],int sersp[],double VBSC[],double NALFn0, double theta,double fai,double Sigmam[][4])
{
	int nmax1=n+5;
	double *NPAI=(double*)malloc(nmax1*sizeof(double));
	double *NTAO=(double*)malloc(nmax1*sizeof(double));
	Norm_AngularFunc_COL(n,theta,NPAI,NTAO);
	for(int i=0;i<3;i++) {for(int j=0;j<4;j++) Sigmam[i][j]=0.0;}
	for(int m=1;m<=serm[n];m++)
	{
		double *BSCnm=&VBSC[sersp[n]+m*8];
		double cosmfai=cos(m*fai),sinmfai=sin(m*fai);
//		Radial
		Sigmam[0][0]+=NPAI[m]/m*((BSCnm[1]+BSCnm[5])*cosmfai-(BSCnm[2]-BSCnm[6])*sinmfai);
		Sigmam[0][1]+=NPAI[m]/m*((BSCnm[2]+BSCnm[6])*cosmfai+(BSCnm[1]-BSCnm[5])*sinmfai);
//		Theta
		Sigmam[1][0]+=NTAO[m]*((BSCnm[1]+BSCnm[5])*cosmfai-(BSCnm[2]-BSCnm[6])*sinmfai);
		Sigmam[1][1]+=NTAO[m]*((BSCnm[2]+BSCnm[6])*cosmfai+(BSCnm[1]-BSCnm[5])*sinmfai);
		Sigmam[1][2]+=NPAI[m]*((BSCnm[3]-BSCnm[7])*cosmfai-(BSCnm[4]+BSCnm[8])*sinmfai);
		Sigmam[1][3]+=NPAI[m]*((BSCnm[4]-BSCnm[8])*cosmfai+(BSCnm[3]+BSCnm[7])*sinmfai);
//		Fai
		Sigmam[2][0]+=NPAI[m]*((BSCnm[1]-BSCnm[5])*cosmfai-(BSCnm[2]+BSCnm[6])*sinmfai);
		Sigmam[2][1]+=NPAI[m]*((BSCnm[2]-BSCnm[6])*cosmfai+(BSCnm[1]+BSCnm[5])*sinmfai);
		Sigmam[2][2]+=NTAO[m]*((BSCnm[3]+BSCnm[7])*cosmfai-(BSCnm[4]-BSCnm[8])*sinmfai);
		Sigmam[2][3]+=NTAO[m]*((BSCnm[4]+BSCnm[8])*cosmfai+(BSCnm[3]-BSCnm[7])*sinmfai);
	}
	double *BSCn0=&VBSC[sersp[n]];
	Sigmam[0][0]=BSCn0[1]*NALFn0+Sigmam[0][0]*sin(theta);
	Sigmam[0][1]=BSCn0[2]*NALFn0+Sigmam[0][1]*sin(theta);
	Sigmam[1][0]+=BSCn0[1]*NTAO[0];
	Sigmam[1][1]+=BSCn0[2]*NTAO[0];
	Sigmam[2][2]+=BSCn0[3]*NTAO[0];
	Sigmam[2][3]+=BSCn0[4]*NTAO[0];
	free(NPAI),free(NTAO); NPAI=NTAO=NULL;
	return;
}


//	***************************************************************************************
//	***************************************************************************************
//	Mie Coefficients
//	***************************************************************************************
//	a[n]=> (CMIE[8*n-7], CMIE[8*n-6])		b[n]=> (CMIE[8*n-5], CMIE[8*n-4])
//	c[n]=> (CMIE[8*n-3], CMIE[8*n-2])		d[n]=> (CMIE[8*n-1], CMIE[8*n-0])
//	***************************************************************************************
//	***************************************************************************************
void Modified_CMIE(int nmax,double rm,double im,double alpha,double CMIE[])
{
	CMIE[0]=nmax;
	if(im<0.0) im=fabs(im);
	int nmax1=nmax+5;
	int n;
	double *D1x=(double*)calloc(nmax1,sizeof(double));
	double *reD1y=(double*)calloc(nmax1,sizeof(double)),*imD1y=(double*)calloc(nmax1,sizeof(double));
	double *reD3x=(double*)calloc(nmax1,sizeof(double)),*imD3x=(double*)calloc(nmax1,sizeof(double));
	double *reR13x=(double*)calloc(nmax1,sizeof(double)),*imR13x=(double*)calloc(nmax1,sizeof(double));
	double *reQ11xy=(double*)calloc(nmax1,sizeof(double)),*imQ11xy=(double*)calloc(nmax1,sizeof(double));

	double comx[2],comy[2],ResLENTZ[2],rma=rm*alpha,ima=im*alpha,nalpha;
//	D1n(x) & D1n(y):****************************************************************
	LENTZ(nmax,alpha,0.0,ResLENTZ); D1x[nmax]=ResLENTZ[0];
	LENTZ(nmax,rma,ima,ResLENTZ); reD1y[nmax]=ResLENTZ[0],imD1y[nmax]=ResLENTZ[1];
	for(n=nmax;n>0;n--)
	{
		nalpha=n/alpha,D1x[n-1]=nalpha-1.0/(nalpha+D1x[n]);
		CplxDiv(n,0.0,rma,ima,comx);
		CplxDiv(1.0,0.0,comx[0]+reD1y[n],comx[1]+imD1y[n],comy);
		reD1y[n-1]=comx[0]-comy[0],imD1y[n-1]=comx[1]-comy[1];
	}
//	D3n(x) & R13n(x):***************************************************************
	reD3x[0]=0.0,reR13x[0]=pow(sin(alpha),2);
	imD3x[0]=1.0,imR13x[0]=0.5*sin(2.0*alpha);
	for(n=1;n<=nmax;n++)
	{
		nalpha=n/alpha;
		CplxDiv(1.0,0.0,nalpha-reD3x[n-1],-imD3x[n-1],comx);
		reD3x[n]=comx[0]-nalpha,imD3x[n]=comx[1];
		CplxMul(reR13x[n-1],imR13x[n-1],nalpha+reD3x[n],imD3x[n],comx);
		reR13x[n]=comx[0]/(nalpha+D1x[n]),imR13x[n]=comx[1]/(nalpha+D1x[n]);
	}
//	Q(x,y):*************************************************************************
	double exp2Nimy=exp(-2*ima);
	CplxDiv(2.0*sin(alpha),0.0,sin(rma)*(1.0+exp2Nimy),cos(rma)*(1.0-exp2Nimy),comx);

	reQ11xy[0]=comx[0],imQ11xy[0]=comx[1];
	for(n=1;n<=nmax;n++)
	{
		double f=n/alpha+D1x[n];
		CplxDiv(n,0.0,rma,ima,comy);
		CplxMul(reQ11xy[n-1],imQ11xy[n-1],comy[0]+reD1y[n],comy[1]+imD1y[n],comx);
		reQ11xy[n]=comx[0]/f,imQ11xy[n]=comx[1]/f;
	}

//	********************************************************************************
	double real_tmp1,imag_tmp1,real_tmp2,imag_tmp2,real_tmp3,imag_tmp3;
	double MIEan[2],MIEbn[2],MIEcn[2],MIEdn[2];
	for(n=1;n<=nmax;n++)
	{
//	D1n[y]/m & m*D1n[y]	************************************************************
		CplxDiv(reD1y[n],imD1y[n],rm,im,comx); real_tmp1=comx[0],imag_tmp1=comx[1];
		CplxMul(reD1y[n],imD1y[n],rm,im,comx); real_tmp2=comx[0],imag_tmp2=comx[1];
		CplxMul(reQ11xy[n],imQ11xy[n],reD3x[n]-D1x[n],imD3x[n],comx); real_tmp3=comx[0],imag_tmp3=comx[1];
//	an	****************************************************************************
		CplxDiv(real_tmp1-D1x[n],imag_tmp1,real_tmp1-reD3x[n],imag_tmp1-imD3x[n],comx);
		CplxMul(reR13x[n],imR13x[n],comx[0],comx[1],MIEan);
//	bn	****************************************************************************
		CplxDiv(real_tmp2-D1x[n],imag_tmp2,real_tmp2-reD3x[n],imag_tmp2-imD3x[n],comx);
		CplxMul(reR13x[n],imR13x[n],comx[0],comx[1],MIEbn);
//	cn	****************************************************************************
		CplxDiv(reD1y[n],imD1y[n],rm,im,comx);
		CplxDiv(real_tmp3,imag_tmp3,reD3x[n]-comx[0],imD3x[n]-comx[1],MIEcn);
//	dn	****************************************************************************
		CplxDiv(reD3x[n],imD3x[n],rm,im,comx);
		CplxDiv(real_tmp3,imag_tmp3,comx[0]-reD1y[n],comx[1]-imD1y[n],MIEdn);
//	********************************************************************************
		CMIE[8*n-7]=MIEan[0],CMIE[8*n-6]=MIEan[1],CMIE[8*n-5]=MIEbn[0],CMIE[8*n-4]=MIEbn[1];
		CMIE[8*n-3]=MIEcn[0],CMIE[8*n-2]=MIEcn[1],CMIE[8*n-1]=MIEdn[0],CMIE[8*n]=MIEdn[1];
	}
	free(reD1y),free(imD1y),free(D1x),free(reD3x),free(imD3x),free(reR13x),free(imR13x),free(reQ11xy),free(imQ11xy);
	reD1y=imD1y=D1x=reD3x=imD3x=reR13x=imR13x=reQ11xy=imQ11xy=NULL;
	return;
}


//	***************************************************************************************
//	***************************************************************************************
//		Read the BSCs:
//		serm[0]=nmax		serm[n]=serm[nr]=mr=mmax
//		for specific n: m=0 -> mmax=serm[n]
//		A(n,+m):	reA[n,+m]=VBSC[sersp[n]+m*8+1]		imA[n,+m]=VBSC[sersp[n]+m*8+2]
//		B(n,+m):	reB[n,+m]=VBSC[sersp[n]+m*8+3]		imB[n,+m]=VBSC[sersp[n]+m*8+4]
//		A(n,-m):	reA[n,-m]=VBSC[sersp[n]+m*8+5]		imA[n,-m]=VBSC[sersp[n]+m*8+6]
//		B(n,-m):	reB[n,-m]=VBSC[sersp[n]+m*8+7]		imB[n,-m]=VBSC[sersp[n]+m*8+8]
//	***************************************************************************************
//	***************************************************************************************
int Norm_BSC_Reading(int ctrl,char Pfile0[],char Pfile1[],int nmin,int nmax,int serm[],int sersp[],double VBSC[])
{
	FILE *fp;
	int n,m,nr,mr,sr,k;
	if(ctrl==1)
	{
		for(n=1;n<=nmin;n++) sersp[n]=serm[n]=0;
		fp=fopen(Pfile1,"r");
		for(n=nmin;n<=nmax;n++)
		{
			fscanf(fp,"%d\t%d\t%d\n",&nr,&mr,&sr);
			if(nr==0 && mr==0 && sr==0) {nmax=n-1;break;}
			if(nr!=n || mr<1 || sr!=sersp[n])
			{
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
				for(int warn=1;warn<=3;warn++) printf("\n\n\tBSC may be wrong, requare recalculation...\n\n"),getche();
				SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
				return 0;
			}
			else
			{
				serm[nr]=mr;
				sersp[n+1]=sersp[n]+(mr+1)*8;
			}
		}
		serm[0]=nmax;
		fclose(fp);
	}
	else
	{
		VBSC[0]=sersp[nmax+1];
		fp=fopen(Pfile0,"r");
		for(k=0,n=nmin;n<=serm[0];n++)
		{
			for(m=0;m<=serm[n];m++)
			{
				fscanf(fp,"%d\t%d",&nr,&mr);
				for(int i=0;i<8;i++) {k++;fscanf(fp,"\t%le",&VBSC[k]);}
				fscanf(fp,"\n");
			}
		}
		fclose(fp);
	}
	return 1;
}



//	***************************************************************************************
//	***************************************************************************************
//	calculate the BSCs and store the result in a file
//	***************************************************************************************
//	***************************************************************************************
int Norm_BSC_Writing(int method,char Pfile0[],char Pfile1[],int nmin,int nmax,double DataBeam[],int serm[])
{
	double APnm[3],BPnm[3],ANnm[3],BNnm[3];
	double APcm[3],BPcm[3],ANcm[3],BNcm[3];
	int RecNumBSC=0;
	char Pfilet[100];
	if(method==0) strcpy(Pfilet,"ASD_CPUT.txt");
	else if(method==1) strcpy(Pfilet,"Quad_CPUT.txt");
	else if(method==3) strcpy(Pfilet,"piASD_CUPT.txt");
	else strcpy(Pfilet,"LAM_CPUT.txt");

	FILE *fp,*fp1,*fp2;
	fp=fopen(Pfile0,"w"),fp1=fopen(Pfile1,"w"),fp2=fopen(Pfilet,"w");
	serm[0]=nmax;
	double ruo0=Funcabs(DataBeam[4],DataBeam[5]);

	for(int n=nmin;n<=nmax;n++)
	{
		clock_t start,finish;
		start=clock();
		DataBeam[0]=n+0.5;
		printf("\n\n****************************************************************************\n");
		int m=0;
		if(method<=1)
		{
			if(ruo0==0.0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=0.0;
			else Norm_BSC_Int(method,n,m,DataBeam,APnm,BPnm);
			printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\tIntegrand\n",nmax,n,m,APnm[2],BPnm[2]);
		}
		else if(method==2)
		{
			if(ruo0==0.0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=0.0;
			else Norm_BSC_LAM(n,m,DataBeam,APnm,BPnm);
			printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\n",nmax,n,m,APnm[2],BPnm[2]);
		}
		else if(method==3)
		{
			if(ruo0==0.0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=0.0;
			else Norm_BSC_Int(method,n,m,DataBeam,APnm,BPnm);
			printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\tIntegrand\n",nmax,n,m,APnm[2],BPnm[2]);
		}
		else
		{
			if(ruo0==0.0)
			{
				APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=0.0;
				APcm[0]=APcm[1]=APcm[2]=BPcm[0]=BPcm[1]=BPcm[2]=0.0;
			}
			else
			{
				Norm_BSC_LAM(n,m,DataBeam,APcm,BPcm);
				Norm_BSC_Int(0,n,m,DataBeam,APnm,BPnm);
			}
			printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\n",nmax,n,m,APnm[2],BPnm[2]);
		}

		if(method<=3) fprintf(fp,"%d\t%d\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n",n,m,APnm[0],APnm[1],BPnm[0],BPnm[1],APnm[0],APnm[1],BPnm[0],BPnm[1]);
		else fprintf(fp,"%d\t%d\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n",n,m,
			APnm[0]-APcm[0],APnm[1]-APcm[1],BPnm[0]-BPcm[0],BPnm[1]-BPcm[1],APnm[0]-APcm[0],APnm[1]-APcm[1],BPnm[0]-BPcm[0],BPnm[1]-BPcm[1]);

		double minratio=1e-5,Vmax=0.0,VCrit,V0=1e-10;
		int NumberFuncCalOld,breakconmold=0,breakconm;

		for(m=1;m<=n;m++)
		{
			if(method<=1)
			{
				NumberFuncCalOld=NumberFuncCal;
				if(ruo0==0.0 && m%2==0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=ANnm[0]=ANnm[1]=ANnm[2]=BNnm[0]=BNnm[1]=BNnm[2]=0.0;
				else
				{
					Norm_BSC_Int(method,n,m,DataBeam,APnm,BPnm);
					Norm_BSC_Int(method,n,-m,DataBeam,ANnm,BNnm);
				}
				printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\t%d\n",nmax,n,m,APnm[2],BPnm[2],NumberFuncCal-NumberFuncCalOld);
			}
			else if(method==3)
			{
				NumberFuncCalOld=NumberFuncCal;
				if(ruo0==0.0 && m%2==0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=ANnm[0]=ANnm[1]=ANnm[2]=BNnm[0]=BNnm[1]=BNnm[2]=0.0;
				else
				{
					Norm_BSC_Int(method,n,m,DataBeam,APnm,BPnm);
					Norm_BSC_Int(method,n,-m,DataBeam,ANnm,BNnm);
				}
				printf("\rnmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\t%d\n",nmax,n,m,APnm[2],BPnm[2],NumberFuncCal-NumberFuncCalOld);
			}
			else if(method==2)
			{
				if(ruo0==0.0 && m%2==0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=ANnm[0]=ANnm[1]=ANnm[2]=BNnm[0]=BNnm[1]=BNnm[2]=0.0;
				else
				{
					Norm_BSC_LAM(n,m,DataBeam,APnm,BPnm);
					Norm_BSC_LAM(n,-m,DataBeam,ANnm,BNnm);
				}
				printf("nmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\n",nmax,n,m,APnm[2],BPnm[2]);
			}
			else
			{
				if(ruo0==0.0 && m%2==0) APnm[0]=APnm[1]=APnm[2]=BPnm[0]=BPnm[1]=BPnm[2]=ANnm[0]=ANnm[1]=ANnm[2]=BNnm[0]=BNnm[1]=BNnm[2]=0.0;
				else
				{
					Norm_BSC_LAM(n,m,DataBeam,APcm,BPcm);
					Norm_BSC_LAM(n,-m,DataBeam,ANcm,BNcm);
					Norm_BSC_Int(0,n,m,DataBeam,APnm,BPnm);
					Norm_BSC_Int(0,n,-m,DataBeam,ANnm,BNnm);
				}
				printf("nmax=%d\tn=%d\tm=%d\t%.3e\t%.3e\n",nmax,n,m,APnm[2],BPnm[2]);
			}
			if(method<=3) fprintf(fp,"%d\t%d\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n",n,m,APnm[0],APnm[1],BPnm[0],BPnm[1],ANnm[0],ANnm[1],BNnm[0],BNnm[1]);
			else fprintf(fp,"%d\t%d\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n",n,m,
				APnm[0]-APcm[0],APnm[1]-APcm[1],BPnm[0]-BPcm[0],BPnm[1]-BPcm[1],ANnm[0]-ANcm[0],ANnm[1]-ANcm[1],BNnm[0]-BNcm[0],BNnm[1]-BNcm[1]);
			if(Vmax<APnm[2]) Vmax=APnm[2],VCrit=minratio*Vmax;
			if(Vmax<BPnm[2]) Vmax=BPnm[2],VCrit=minratio*Vmax;
			if(Vmax<ANnm[2]) Vmax=ANnm[2],VCrit=minratio*Vmax;
			if(Vmax<BNnm[2]) Vmax=BNnm[2],VCrit=minratio*Vmax;
			if(VCrit<V0) VCrit=V0;
//			if ((APnm[2]<VCrit && BPnm[2]<VCrit && ANnm[2]<VCrit && BNnm[2]<VCrit) || Vmax==0.0) {serm[n]=m;break;}
			breakconm=((APnm[2]<VCrit && BPnm[2]<VCrit && ANnm[2]<VCrit && BNnm[2]<VCrit) || Vmax==0.0)? 1:0;
			if(breakconmold==1 && breakconm==1) {serm[n]=m;break;}
			else breakconmold=breakconm;
		}
		finish=clock();
		double Total_time=(double)(finish-start)/CLOCKS_PER_SEC;
		fprintf(fp2,"%d\t%d\t%le\t%d\n",n,m,Total_time,NumberFuncCal);
		if(m<n) {fprintf(fp1,"%d\t%d\t%d\n",n,m,RecNumBSC);RecNumBSC+=(m+1)*8;}
		else {fprintf(fp1,"%d\t%d\t%d\n",n,n,RecNumBSC);serm[n]=n,RecNumBSC+=(n+1)*8;}
		fclose(fp),fclose(fp1),fclose(fp2);
		fp=fopen(Pfile0,"a"),fp1=fopen(Pfile1,"a"),fp2=fopen(Pfilet,"a");
	}
	fclose(fp),fclose(fp1),fclose(fp2);
	return 0;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	BSC calculation with the Localized Approximation
//	**************************************************************************************************
//	**************************************************************************************************
void Norm_BSC_LAM(int n,int m,double DataBeam[],double Anm[],double Bnm[])
{
	int absm=abs(m);
	if(n<absm || n==0) {Anm[0]=Anm[1]=Anm[2]=Bnm[0]=Bnm[1]=Bnm[2]=0.0;return;}
//	*******************************************************************************
	double ruon=n+0.5,sx=DataBeam[2],sy=DataBeam[3],X0=DataBeam[4],Y0=DataBeam[5],Z0=DataBeam[6];
	double sx2=sx*sx,sy2=sy*sy,sx2X0=sx2*X0,sy2Y0=sy2*Y0;
//	*******************************************************************************
//	Qx & Qy��
	double Qx[2],Qy[2];
	CplxDiv(1.0,0.0,-2.0*sx2*Z0,-1.0,Qx);
	CplxDiv(1.0,0.0,-2.0*sy2*Z0,-1.0,Qy);

//	sqrt(Qx*Qy)��
	double tmp1[2],tmp2[2],tmp3[2],SQxy[2];
	CplxMul(Qx[0],Qx[1],Qy[0],Qy[1],tmp3);
	SQxy[1]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])-tmp3[0]));
	SQxy[0]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])+tmp3[0]));
	if(tmp3[1]<0.0) SQxy[0]*=-1.0;

//	*******************************************************************************
//	Gxy
	double Gxy[2];
	CplxMul(Qx[0]*sx2X0,Qx[1]*sx2X0,Qx[0]*sx2X0,Qx[1]*sx2X0,tmp1);
	CplxMul(Qy[0]*sy2Y0,Qy[1]*sy2Y0,Qy[0]*sy2Y0,Qy[1]*sy2Y0,tmp2);
	tmp3[0]=tmp1[0]+tmp2[0],tmp3[1]=tmp1[1]+tmp2[1];
	Gxy[1]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])-tmp3[0]));
	Gxy[0]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])+tmp3[0]));
	if(tmp3[1]<0.0) Gxy[0]*=-1.0;

//	*******************************************************************************
//	cosKsi & sinksi
	double cosKsi[2],sinKsi[2];
	CplxDiv(Qx[0]*sx2X0,Qx[1]*sx2X0,Gxy[0],Gxy[1],cosKsi);
	CplxDiv(Qy[0]*sy2Y0,Qy[1]*sy2Y0,Gxy[0],Gxy[1],sinKsi);
//	*******************************************************************************
//	2A
	double DblA[2],tmppara=0.5*ruon*ruon;
	DblA[0]=tmppara*(sy2*Qy[1]-sx2*Qx[1]);
	DblA[1]=tmppara*(sx2*Qx[0]-sy2*Qy[0]);
//	*******************************************************************************
//	Beta
	double Beta[2];
	tmppara=2.0*ruon;
	Beta[0]=tmppara*Gxy[1],Beta[1]=-tmppara*Gxy[0];
//	*******************************************************************************
//	fai00
	double fai00[2];
	tmppara=0.5*ruon*ruon;
	tmp1[0]=(Qx[0]*sx2+Qy[0]*sy2)*tmppara;
	tmp1[1]=(Qx[1]*sx2+Qy[1]*sy2)*tmppara;
	tmppara=pow(sx*X0,2);
	tmp1[0]+=Qx[0]*tmppara,tmp1[1]+=Qx[1]*tmppara;
	tmppara=pow(sy*Y0,2);
	tmp1[0]+=Qy[0]*tmppara-Z0,tmp1[1]+=Qy[1]*tmppara;
	tmp2[0]=-tmp1[1],tmp2[1]=tmp1[0];

	if(Beta[0]>=0.0) tmp2[0]+=Beta[0],tmp2[1]+=Beta[1];
	else tmp2[0]-=Beta[0],tmp2[1]-=Beta[1];
	if(DblA[0]>=0.0 ) tmp2[0]+=DblA[0],tmp2[1]+=DblA[1];
	else tmp2[0]-=DblA[0],tmp2[1]-=DblA[1];

	tmppara=exp(tmp2[0]);
	CplxMul(SQxy[0],SQxy[1],cos(tmp2[1]),sin(tmp2[1]),tmp1);
	fai00[0]=tmppara*tmp1[1],fai00[1]=-tmppara*tmp1[0];

//	*******************************************************************************
	double FImAnm[2],FImBnm[2];
	FuncRelaModifiedBesselLAM(m,DblA,Beta,cosKsi,sinKsi,FImAnm,FImBnm);
//	*******************************************************************************
	CplxMul(fai00[0],fai00[1],FImAnm[0],FImAnm[1],tmp1);
	CplxMul(fai00[0],fai00[1],FImBnm[0],FImBnm[1],tmp2);
	if(m==0)
	{
		double eznm=-n*(n+1.0)/pow(n+0.5,1.5);
		CplxMul(0.0,eznm,tmp1[0],tmp1[1],Anm);
		CplxMul(-eznm,0.0,tmp2[0],tmp2[1],Bnm);
	}
	else
	{
		double paraTZ[2],znm=sqrt(n+0.5);
		for(int k=1;k<=absm;k++) znm*=sqrt((n+k)*(n-k+1.0))/(n+0.5);
		switch((absm-1)%4)
		{
			case 0: paraTZ[0]=-znm,paraTZ[1]=0.0; break;
			case 1: paraTZ[0]=0.0,paraTZ[1]=znm; break;
			case 2: paraTZ[0]=znm,paraTZ[1]=0.0; break;
			case 3: paraTZ[0]=0.0,paraTZ[1]=-znm; break;
		}
		CplxMul(paraTZ[0],paraTZ[1],tmp1[0],tmp1[1],Anm);
		CplxMul(paraTZ[0],paraTZ[1],-tmp2[1],tmp2[0],Bnm);
	}
	Anm[2]=Funcabs(Anm[0],Anm[1]),Bnm[2]=Funcabs(Bnm[0],Bnm[1]);
	return;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	For BSC LAM
//	**************************************************************************************************
//	**************************************************************************************************
void FuncRelaModifiedBesselLAM(int m,double DblA[],double Beta[],double cosKsi[],double sinKsi[],double FImAnm[2],double FImBnm[2])
{
	double Vmax=0.0,absDblA=Funcabs(DblA[0],DblA[1]);
	double absBeta=Funcabs(Beta[0],Beta[1]);
	if(absDblA==0.0)	//	circular Gaussian beam
	{
		double resmn1[2],resm0[2],resmp1[2],tmp1[2],tmp2[2];
		CMPLE_EqBesselI_3Order(m,Beta[0],Beta[1],resmn1,resm0,resmp1);

		int mn=m-1;
		if(mn>=0) CplxPow(mn,cosKsi[0]+sinKsi[1],cosKsi[1]-sinKsi[0],tmp1);
		else CplxPow(-mn,cosKsi[0]-sinKsi[1],cosKsi[1]+sinKsi[0],tmp1);
		CplxMul(resmn1[0],resmn1[1],tmp1[0],tmp1[1],tmp1);
		double Vtemp1=Funcabs(tmp1[0],tmp1[1]);

		mn=m+1;
		if(mn>=0) CplxPow(mn,cosKsi[0]+sinKsi[1],cosKsi[1]-sinKsi[0],tmp2);
		else CplxPow(-mn,cosKsi[0]-sinKsi[1],cosKsi[1]+sinKsi[0],tmp2);
		CplxMul(resmp1[0],resmp1[1],tmp2[0],tmp2[1],tmp2);
		double Vtemp2=Funcabs(tmp2[0],tmp2[1]);
		Vmax=(Vtemp1>Vtemp2)?Vtemp1:Vtemp2;

		FImAnm[0]=tmp2[0]+tmp1[0],FImAnm[1]=tmp2[1]+tmp1[1];
		FImBnm[0]=tmp2[0]-tmp1[0],FImBnm[1]=tmp2[1]-tmp1[1];
	}
	else				//	EGB
	{
		double tmp1[2],FImAnmtmp[2],FImBnmtmp[2];
		int p,s,absp,abss,abssp1,smax,pmax,s1,s2,p1,p2,absm=abs(m);
		pmax=int(absBeta+10.5*pow(absBeta,0.341)+20.5);
		s1=-(pmax+m)/2,s2=(pmax-m)/2;
		if(m>0 && s2<20) s2=20;
		if(m<0 && s1>-20) s1=-20;
		smax=((s1+s2)<0)?-s1:s2;
		if((p1=m+2*s1)<0) p1=-p1;
		if((p2=m+2*s2)<0) p2=-p2;
		pmax=(p1>p2)?p1:p2;

		double *reEqBslI2A=(double*)calloc((smax+5),sizeof(double));
		double *imEqBslI2A=(double*)calloc((smax+5),sizeof(double));
		double *reEqBslIBeta=(double*)calloc((pmax+5),sizeof(double));
		double *imEqBslIBeta=(double*)calloc((pmax+5),sizeof(double));
		CMPLE_EqBesselI(smax+3,DblA[0],DblA[1],reEqBslI2A,imEqBslI2A);
		CMPLE_EqBesselI(pmax+3,Beta[0],Beta[1],reEqBslIBeta,imEqBslIBeta);

//	=====================================================================================================
		FImAnm[0]=FImAnm[1]=0.0;
		FImBnm[0]=FImBnm[1]=0.0;
		double Vtmp,Vmax=0.0;
		for(s=s1;s<=s2;s++)
		{
			abss=abs(s);
			p=s+1,abssp1=abs(p);
			p=m+2*s+1,absp=abs(p);
			if(p>=0) CplxPow(p,cosKsi[0]+sinKsi[1],cosKsi[1]-sinKsi[0],tmp1);
			else CplxPow(-p,cosKsi[0]-sinKsi[1],cosKsi[1]+sinKsi[0],tmp1);
			CplxMul(reEqBslIBeta[absp],imEqBslIBeta[absp],tmp1[0],tmp1[1],tmp1);

			CplxMul(reEqBslI2A[abss]+reEqBslI2A[abssp1],imEqBslI2A[abss]+imEqBslI2A[abssp1],tmp1[0],tmp1[1],FImAnmtmp);
			CplxMul(reEqBslI2A[abss]-reEqBslI2A[abssp1],imEqBslI2A[abss]-imEqBslI2A[abssp1],tmp1[0],tmp1[1],FImBnmtmp);
			FImAnm[0]+=FImAnmtmp[0],FImAnm[1]+=FImAnmtmp[1];
			FImBnm[0]+=FImBnmtmp[0],FImBnm[1]+=FImBnmtmp[1];
			if(Vmax<(Vtmp=Funcabs(FImAnmtmp[0],FImAnmtmp[1]))) Vmax=Vtmp;
			if(Vmax<(Vtmp=Funcabs(FImBnmtmp[0],FImBnmtmp[1]))) Vmax=Vtmp;
		}
		free(reEqBslI2A),free(imEqBslI2A),free(reEqBslIBeta),free(imEqBslIBeta);
		reEqBslI2A=imEqBslI2A=reEqBslIBeta=imEqBslIBeta=NULL;
		if(Funcabs(FImAnm[0],FImAnm[1])<1e-14*Vmax) FImAnm[0]=FImAnm[1]=0.0;
		if(Funcabs(FImBnm[0],FImBnm[1])<1e-14*Vmax) FImBnm[0]=FImBnm[1]=0.0;
	}
	return;
}





//	**************************************************************************************************
//	**************************************************************************************************
//	BSC calculation with Quadreture or ASD
//	**************************************************************************************************
//	**************************************************************************************************
void Norm_BSC_Int(int intway,int n,int m,double DataBeam[],double Anm[],double Bnm[])
{
	int absm=abs(m);
	if(n<absm || n==0)
	{
		Anm[0]=Anm[1]=Anm[2]=Bnm[0]=Bnm[1]=Bnm[2]=0.0;
		return;
	}
//	*******************************************************************************
	double Q0[2];
	if(intway==1)	//	Quadrature method
	{
		double ruon=DataBeam[0]=n+0.5,Z0=DataBeam[6];
		double RB0=cos(ruon),RB1=sin(ruon),RB2;
		for(int i=1;i<=n;i++) RB2=(2*i-1)*RB1/ruon-RB0,RB0=RB1,RB1=RB2;
		double paraQ0=-ruon*ruon/(2.0*n+1.0)/RB2;

		switch(n%4)
		{
			case 0: Q0[0]=paraQ0,Q0[1]=0.0;break;
			case 1: Q0[0]=0.0,Q0[1]=-paraQ0; break;
			case 2: Q0[0]=-paraQ0,Q0[1]=0.0; break;
			case 3: Q0[0]=0.0,Q0[1]=paraQ0; break;
		}
		CplxMul(Q0[0],Q0[1],cos(Z0),-sin(Z0),Q0);
		double ABnm[4];
		Integral_Q(intway,n,m,DataBeam,ABnm);
		CplxMul(Q0[0],Q0[1],ABnm[0],ABnm[1],Anm);
		CplxMul(Q0[0],Q0[1],ABnm[2],ABnm[3],Bnm);
	}
//	*******************************************************************************
	else if(intway==0)		//	angular spectra decomposition
	{
		double ABnm[4],tmp[2];
		int absm=abs(m);
		double sx=DataBeam[2],sy=DataBeam[3];
		double paraQ0=-1.0/(2.0*n+1.0)/sx/sy;
		switch((absm*4+1-m)%4)
		{
			case 0: tmp[0]=paraQ0,tmp[1]=0.0;break;
			case 1: tmp[0]=0.0,tmp[1]=paraQ0; break;
			case 2: tmp[0]=-paraQ0,tmp[1]=0.0; break;
			case 3: tmp[0]=0.0,tmp[1]=-paraQ0; break;
		}
		Integral_Q(intway,n,m,DataBeam,ABnm);
		CplxMul(tmp[0],tmp[1],ABnm[0],ABnm[1],Anm);
		CplxMul(tmp[1],-tmp[0],ABnm[2],ABnm[3],Bnm);
	}
//	*******************************************************************************
	else	//	piASD-BSCs
	{
		_Complex double fact_tm;
		_Complex double Unm;
		double res1[2],res2[2];
		_Complex double tm1,tm2,Tm1,Tm2;
		int l_order=DataBeam[9],p_order=DataBeam[10];
		double k=DataBeam[1],s=DataBeam[2],z0=DataBeam[16],x0=DataBeam[14],y0=DataBeam[15];
		double Rloc=n+0.5;
		double a1=2.0*pow(s,2)*k*z0;
		double kosai0=atan2(y0,x0);
		fact_tm=(-1/(2*pi))*cpow(1+I*a1,p_order)/cpow(1-I*a1,(l_order+p_order+1));	//fact of Tm
		
		Unm=cpow(I,1-m)*pow(Rloc,0.5-abs(m))*sqrt(def(n+abs(m),n-abs(m)));
		Integraltm(n,m-1,DataBeam,res1);
		Integraltm(n,m+1,DataBeam,res2);
		tm1=res1[0]+I*res1[1];
		tm2=res2[0]+I*res2[1];
		Tm1=fact_tm*tm1;
		Tm2=fact_tm*tm2;
		_Complex double A=Unm*cexp(-1*I*k*z0)*(cexp(-1*I*(m-l_order-1)*kosai0)*Tm1+cexp(-1*I*(m-l_order+1)*kosai0)*Tm2);
		_Complex double B=Unm*cexp(-1*I*k*z0)*(cexp(-1*I*(m-l_order-1)*kosai0)*Tm1-cexp(-1*I*(m-l_order+1)*kosai0)*Tm2)/I;
		Anm[0]=creal(A);Anm[1]=cimag(A);
		Bnm[0]=creal(B);Bnm[1]=cimag(B);
	}

	Anm[2]=Funcabs(Anm[0],Anm[1]),Bnm[2]=Funcabs(Bnm[0],Bnm[1]);
	return;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Integration:tm
//	**************************************************************************************************
//	**************************************************************************************************
void Integraltm(int n,int m,double DataBeam[],double res[]){
	double s=DataBeam[2],X0=DataBeam[4],Y0=DataBeam[5],wavenumber0=DataBeam[1],z0=DataBeam[16];
	int l_order=DataBeam[9],p_order=DataBeam[10];
	double r0,Rloc,a1;
	r0=sqrt(X0*X0+Y0*Y0);
	Rloc=n+0.5;
	a1=2*pow(s,2)*wavenumber0*z0;
	
	int kmax=80,k,i,j;
	double lowlimit=0.0,uplimit=2.0*pi,length=uplimit-lowlimit,tmp[3];
	res[0]=res[1]=0.0;
	double step=(uplimit-lowlimit)/kmax;
	double IntFuncmax=0.0;
//	printf("%f\n",step);
	for(k=1;k<=kmax;k++){
		double a=lowlimit+(k-1)*step,b=lowlimit+k*step;
		int Int_ok=Romberg(n,m,a,b,DataBeam,tmp);
		for(j=0;j<=1;j++)
			res[j]+=tmp[j];
		if(IntFuncmax<tmp[2])
			IntFuncmax=tmp[2];
		if(IntFuncmax*1e-10>tmp[2])
			break;
	}
//	printf("%f\t%f\t",res[0],res[1]);
}

// Romberg����
int Romberg(int n,int m,double a,double b,double DataBeam[],double res[]){
	double s=DataBeam[2],X0=DataBeam[4],Y0=DataBeam[5],wavenumber0=DataBeam[1],z0=DataBeam[16];
	int l_order=DataBeam[9],p_order=DataBeam[10];
	double r0,Rloc,a1;
	r0=sqrt(X0*X0+Y0*Y0);
	Rloc=n+0.5;
	a1=2*pow(s,2)*wavenumber0*z0;

	double h=b-a;					// ��������
	const double eps=1e-5;
	const int maxk=10;
	int i,j,k,ep;
	double res1[3],res2[3],R[maxk+2][2],RZ[2],Rp[2],cum[2];
	
//	fun(a,res1);
//	fun(b,res1);
	Integraltm0(a,n,m,DataBeam,res1);
	Integraltm0(b,n,m,DataBeam,res2);
	
	res[2]=(res1[2]<res2[2])? res2[2]:res1[2]; 
	
	for(i=0;i<=1;i++)
		R[0][i]=h*0.5*(res1[i]+res2[i]);
	double tmold,tmnew=Funcabs(R[0][0],R[0][1]);
	
	for(k=1;k<=maxk;k++){
		h*=0.5,ep=1<<(k-1);
		for(i=0;i<=1;i++)
			cum[i]=0.0;
			
		for(j=1;j<=ep;j++){
//			fun(a+(2*j-1)*h,res1);
			Integraltm0(a+(2*j-1)*h,n,m,DataBeam,res1);
			if(res[2]<res1[2])
				res[2]=res1[2];
			for(i=0;i<=1;i++)
				cum[i]+=res1[i];
		}
		for(i=0;i<=1;i++){
			RZ[i]=R[0][i];
			R[0][i]=h*cum[i]+0.5*RZ[i];
		}
		for(j=1;j<=k;j++){
			for(i=0;i<=1;i++){
				Rp[i]=R[j-1][i]+(R[j-1][i]-RZ[i])/(pow(4,j)-1);
				RZ[i]=R[j][i];
				R[j][i]=Rp[i];
			}
		}
		tmold=tmnew;
		tmnew=Funcabs(R[k][0],R[k][1]);
		if(k>4 && (fabs(tmold-tmnew)<=eps*tmnew || tmnew<1e-300)){
		  	for(i=0;i<=1;i++)
				res[i]=R[k][i];
//		  	printf("%f\t%f\n",res[0],res[1]);
			return 1;
		}
	}
	for(i=0;i<=1;i++)
		res[i]=R[maxk][i];
//	printf("%f\t%f\t%f\t%f\t%f\t%f\n",res1[0],res1[1],res1[2],res2[0],res2[1],res2[2]);
//	printf("%f\t%f\n",res[0],res[1]);
	return 0;
}

// integrand �������� 
void Integraltm0(double gamma,int n,int m,double DataBeam[],double res[]){
	double s=DataBeam[2],X0=DataBeam[4],Y0=DataBeam[5],wavenumber0=DataBeam[1],z0=DataBeam[16];
	int l_order=DataBeam[9],p_order=DataBeam[10];
	double r0,Rloc,a1;
	r0=sqrt(X0*X0+Y0*Y0);
	Rloc=n+0.5;
	a1=2*s*s*wavenumber0*z0;
	
	double w1=sqrt(pow(Rloc,2)+pow(r0,2)-2*Rloc*r0*cos(gamma));
	double kappa=asin(r0*sin(gamma)/w1);
	double fa=2*pow(w1*s,2)/(1+pow(a1,2));
	_Complex double fun=cexp(I*l_order*kappa+I*(l_order-m)*gamma-(pow(w1*s,2)/(1-I*a1)))*pow(sqrt(2)*w1*s,l_order)*Laguerre(l_order,p_order,fa);
	res[0]=creal(fun);
	res[1]=cimag(fun);
	res[2]=Funcabs(res[0],res[1]);
}
// ����Unm ����
double def(int fx,int fy){
	double fact=1.0;
	int i;
	if(fx==fy)
		return 1.0;
	else if(fx>fy){
		for(i=1;i<=fx-fy;i++)
			fact=fact*(fy+i);
		return fact;
	}
	else if(fx<fy){
		for(i=1;i<=fy-fx;i++)
			fact=fact*(fx+i);
		return 1.0/fact;
	}
}
//�׳� 
int factorial(double x){
    int i;
    int result=1;
    if(x<0)
        result=0;
    else
        for(i=1; i<=x; ++i)
            result *= i;
    return result;
}

//Laguerre����ʽ 
double Laguerre(double l_order,double p_order,double x){
    double tmp=0.0;
    int m;
    for(m=0;m<=p_order;m++)
        tmp=tmp+pow(-1.0,m)*factorial(l_order+p_order)*pow(x,m)/(factorial(p_order-m)*factorial(l_order+m)*factorial(m));
    return tmp;
}



//	**************************************************************************************************
//	**************************************************************************************************
//	Integration:
//	**************************************************************************************************
//	**************************************************************************************************
void Integral_Q(int intway,int n,int m,double DataBeam[],double res[])
{
	res[0]=res[1]=res[2]=res[3]=0.0;
	double lowlimit=0.0,uplimit=0.5*pi,length=uplimit-lowlimit,tmp[6];
	int absm=abs(m),kmax=(n<20)?20:(n+1);
	int kcr=int(kmax*m/(n+1.0));
	double step=0.5*pi/kmax;
	double IntFuncmax=0.0;
//	**************************************************************************************************
	if(intway==0)	//	ASD
	{
		for(int k=1;k<=kmax;k++)
		{
			double a=lowlimit+(k-1)*step,b=lowlimit+k*step;
			int Count=NumberFuncCal;
			int Int_ok=RombergTrapezoi_SubQ(intway,n,m,a,b,DataBeam,tmp);
			for(int j=0;j<4;j++) res[j]+=tmp[j];
			printf("\rInt:\t%d\t%d\t%.3e\t%.3e\t%d\t",k,kmax,a,b,NumberFuncCal-Count);
			if(IntFuncmax<tmp[4]) IntFuncmax=tmp[4];
			if(IntFuncmax*1e-10>tmp[4]) break;
		}
	}
//	**************************************************************************************************
	else	//	Quadrature
	{
		for(int k=kmax;k>=1;k--)
		{
			int j;
			int Count=NumberFuncCal;
//		0 to 90 degrees
			double a=lowlimit+(k-1)*step,b=lowlimit+k*step;
			int Int_ok=RombergTrapezoi_SubQ(intway,n,m,a,b,DataBeam,tmp);
			for(j=0;j<4;j++) res[j]+=tmp[j];
			if(IntFuncmax<tmp[4]) IntFuncmax=tmp[4];
//		90 to 180 degrees
			double c=pi-b,d=pi-a;
			Int_ok=RombergTrapezoi_SubQ(intway,n,m,c,d,DataBeam,tmp);
			for(j=0;j<4;j++) res[j]+=tmp[j];
			printf("\rint:\t%d\t%d\t%.3e\t%.3e\t%d\t",k,kmax,a,b,NumberFuncCal-Count);
			if(IntFuncmax<tmp[4]) IntFuncmax=tmp[4];
			if(IntFuncmax*1e-10>tmp[4] && k<kcr) break;
		}
	}
	return;
}




//	**************************************************************************************************
//	**************************************************************************************************
//	Romberg integration
//	**************************************************************************************************
//	**************************************************************************************************
int RombergTrapezoi_SubQ(int intway,int n,int m,double a,double b,double DataBeam[],double res[])
{
	res[4]=0.0;
	const double eps=1e-5;
	const int maxk=10;
	int i,j,k,ep,Bessel_ok;
	double res1[6],res2[6],R[maxk+2][4],RZ[4],Rp[4],cum[4],h=b-a;
	if(intway==0)//ASD 
	{
		Integrand_Q0(n,m,DataBeam,a,res1);
		Integrand_Q0(n,m,DataBeam,b,res2);
	}
	else
	{
		Bessel_ok=Integrand_Q1(n,m,DataBeam,a,res1);
		Bessel_ok=Integrand_Q1(n,m,DataBeam,b,res2);
	}
	res[4]=(res1[4]<res2[4])? res2[4]:res1[4];


	for(i=0;i<=3;i++) R[0][i]=h*0.5*(res1[i]+res2[i]);
	double ABnmold,ABnmnew=Funcabs(R[0][0],R[0][1])+Funcabs(R[0][2],R[0][3]);

	for(k=1;k<=maxk;k++)
	{
		h*=0.5,ep=1<<(k-1);
		for(i=0;i<=3;i++) cum[i]=0.0;
		for(j=1;j<=ep;j++)
		{
			if(intway==0) Integrand_Q0(n,m,DataBeam,a+(2*j-1)*h,res1);
			else Bessel_ok=Integrand_Q1(n,m,DataBeam,a+(2*j-1)*h,res1);
			if(res[4]<res1[4]) res[4]=res1[4];
			for(i=0;i<=3;i++) cum[i]+=res1[i];
		}
		for(i=0;i<=3;i++) RZ[i]=R[0][i],R[0][i]=h*cum[i]+0.5*RZ[i];
		for(j=1;j<=k;j++) {for(i=0;i<=3;i++) Rp[i]=R[j-1][i]+(R[j-1][i]-RZ[i])/(pow(4,j)-1),RZ[i]=R[j][i],R[j][i]=Rp[i];}
		ABnmold=ABnmnew,ABnmnew=Funcabs(R[k][0],R[k][1])+Funcabs(R[k][2],R[k][3]);
		if(k>4 && (fabs(ABnmold-ABnmnew)<=eps*ABnmnew || ABnmnew<1e-300)) {for(i=0;i<=3;i++) res[i]=R[k][i]; return 1;}
	}
	for(i=0;i<=3;i++) res[i]=R[maxk][i];
	return 0;
}




//	**************************************************************************************************
//	**************************************************************************************************
//	Integrand for Quadrature method
//	**************************************************************************************************
//	**************************************************************************************************
int Integrand_Q1(int n,int m,double DataBeam[],double theta,double res[])
{
	if(NumberFuncCal<2147483647) NumberFuncCal++;
	else NumberFuncCal=0,NumberHighDigit++;
//	*******************************************************************************
	double ruon=DataBeam[0]=n+0.5,sx=DataBeam[2],sy=DataBeam[3],X0=DataBeam[4],Y0=DataBeam[5],Z0=DataBeam[6];
//	*******************************************************************************
	double sx2=sx*sx,sy2=sy*sy,sx2X0=sx2*X0,sy2Y0=sy2*Y0,sx2X02=sx2X0*sx2X0,sy2Y02=sy2Y0*sy2Y0;
//	*******************************************************************************
	double sintheta=sin(theta),costheta=cos(theta);
	if(sintheta==0.0) {res[0]=res[1]=res[2]=res[3]=0.0;return 1;}
	int absm=abs(m);
	double Pnm=Norm_AssLegFunc_Single(n,absm,theta);
	if(Pnm==0.0) {res[0]=res[1]=res[2]=res[3]=0.0;return 1;}
//	*******************************************************************************
//	Qx & Qy
	double Qx[2],Qy[2];
	double tmppara=2.0*(ruon*costheta-Z0);
	CplxDiv(1.0,0.0,sx2*tmppara,-1.0,Qx);
	CplxDiv(1.0,0.0,sy2*tmppara,-1.0,Qy);
//	*******************************************************************************
//	sqrt(Qx*Qy)
	double tmp1[2],tmp2[2],tmp3[2],SQxy[2];
	CplxMul(Qx[0],Qx[1],Qy[0],Qy[1],tmp3);
	SQxy[1]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])-tmp3[0]));
	SQxy[0]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])+tmp3[0]));
	if(tmp3[1]<0.0) SQxy[0]*=-1.0;

//	*******************************************************************************
//	Gxy
	double Gxy[2];
	CplxMul(Qx[0],Qx[1],Qx[0],Qx[1],tmp1);
	CplxMul(Qy[0],Qy[1],Qy[0],Qy[1],tmp2);
	tmp3[0]=tmp1[0]*sx2X02+tmp2[0]*sy2Y02,tmp3[1]=tmp1[1]*sx2X02+tmp2[1]*sy2Y02;
	Gxy[1]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])-tmp3[0]));
	Gxy[0]=sqrt(0.5*(sqrt(tmp3[0]*tmp3[0]+tmp3[1]*tmp3[1])+tmp3[0]));
	if(tmp3[1]<0.0) Gxy[0]*=-1.0;

//	*******************************************************************************
//	cosKsi & sinksi
	double cosKsi[2],sinKsi[2];
	CplxDiv(Qx[0]*sx2X0,Qx[1]*sx2X0,Gxy[0],Gxy[1],cosKsi);
	CplxDiv(Qy[0]*sy2Y0,Qy[1]*sy2Y0,Gxy[0],Gxy[1],sinKsi);
//	*******************************************************************************
//	2A
	tmppara=0.5*pow(ruon*sintheta,2);
	double DblA[2]={tmppara*(sy2*Qy[1]-sx2*Qx[1]),tmppara*(sx2*Qx[0]-sy2*Qy[0])};
//	*******************************************************************************
//	Beta
	tmppara=2.0*ruon*sintheta;
	double Beta[2]={tmppara*Gxy[1],-tmppara*Gxy[0]};
//	*******************************************************************************
	double FIm[3][2]={0,0,0,0,0,0};
	int Bessel_ok=FuncRelaModifiedBessel(m,DblA,Beta,cosKsi,sinKsi,FIm);
//	*******************************************************************************
	double gamma1[2],gamma2x[2],gamma2y[2],gamma3x[2],gamma3y[2];
//	gamma1
	tmppara=0.5*pow(ruon*sintheta,2);
	tmp1[0]=(Qx[0]*sx2+Qy[0]*sy2)*tmppara;
	tmp1[1]=(Qx[1]*sx2+Qy[1]*sy2)*tmppara;
	tmppara=pow(sx*X0,2);
	tmp1[0]+=Qx[0]*tmppara,tmp1[1]+=Qx[1]*tmppara;
	tmppara=pow(sy*Y0,2);
	tmp1[0]+=Qy[0]*tmppara+ruon*costheta,tmp1[1]+=Qy[1]*tmppara;
	tmp2[0]=-tmp1[1],tmp2[1]=tmp1[0];

	if(Beta[0]>=0.0) tmp2[0]+=Beta[0],tmp2[1]+=Beta[1];
	else tmp2[0]-=Beta[0],tmp2[1]-=Beta[1];
	if(DblA[0]>=0.0) tmp2[0]+=DblA[0],tmp2[1]+=DblA[1];
	else tmp2[0]-=DblA[0],tmp2[1]-=DblA[1];

	tmppara=exp(tmp2[0]);
	CplxMul(SQxy[0],SQxy[1],tmppara*cos(tmp2[1]),tmppara*sin(tmp2[1]),tmp1);
	tmppara=sintheta*Pnm;
	gamma1[0]=tmp1[0]*tmppara,gamma1[1]=tmp1[1]*tmppara;

//	gamma2
	tmppara=2.0*sx2*ruon*costheta;
	gamma2x[0]=sintheta*(1.0-Qx[0]*tmppara);
	gamma2x[1]=-sintheta*Qx[1]*tmppara;
	tmppara=2.0*sy2*ruon*costheta;
	gamma2y[0]=sintheta*(1.0-Qy[0]*tmppara);
	gamma2y[1]=-sintheta*Qy[1]*tmppara;

//	gamma3
	tmppara=4.0*sx2X0*costheta;
	gamma3x[0]=Qx[0]*tmppara,gamma3x[1]=Qx[1]*tmppara;
	tmppara=4.0*sy2Y0*costheta;
	gamma3y[0]=Qy[0]*tmppara,gamma3y[1]=Qy[1]*tmppara;
//	*******************************************************************************
	CplxMul(gamma2x[0],gamma2x[1],FIm[1][0],FIm[1][1],tmp1);
	CplxMul(gamma3x[0],gamma3x[1],FIm[0][0],FIm[0][1],tmp2);
	CplxMul(gamma1[0],gamma1[1],tmp1[0]+tmp2[0],tmp1[1]+tmp2[1],tmp3);
	res[0]=tmp3[0],res[1]=tmp3[1];

	CplxMul(gamma2y[0],gamma2y[1],FIm[2][0],FIm[2][1],tmp1);
	CplxMul(gamma3y[0],gamma3y[1],FIm[0][0],FIm[0][1],tmp2);
	CplxMul(gamma1[0],gamma1[1],tmp2[0]-tmp1[1],tmp2[1]+tmp1[0],tmp3);
	res[2]=tmp3[0],res[3]=tmp3[1];
	for(int i=0;i<4;i++) if(fabs(res[i])<1e-300) res[i]=0.0;
	res[4]=Funcabs(res[0],res[1]),res[5]=Funcabs(res[2],res[3]);
	return 1;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Integrand for ASD
//	**************************************************************************************************
//	**************************************************************************************************
void Integrand_Q0(int n,int m,double DataBeam[],double alpha,double res[])
{
	if(NumberFuncCal<2147483647) NumberFuncCal++;
	else NumberFuncCal=0,NumberHighDigit++;

//	*******************************************************************************
	double sx=DataBeam[2],sy=DataBeam[3],X0=DataBeam[4],Y0=DataBeam[5],Z0=DataBeam[6];
	double sinalpha=sin(alpha),cosalpha=cos(alpha);
	double R0=sqrt(X0*X0+Y0*Y0);
	double Ksi=(R0>0.0)? acos(X0/R0):0.0;
	if(Y0<0.0) Ksi=-Ksi;
	double DblA=(1.0/sy/sy-1.0/sx/sx)*sinalpha*sinalpha/8.0;
	double Beta=R0*sinalpha;
//	*******************************************************************************
	double NPAITAO[2];
	int absm=abs(m);
	Norm_AngularFunc_Single(n,absm,alpha,NPAITAO);
	if(m<0) NPAITAO[0]*=-1;
//	*******************************************************************************
	double Xmn1[3],Xmp1[3],tmp0[2],tmpA[2],tmpB[2];
	FuncXmBesselJB1(m,DblA,Beta,Ksi,Xmn1,Xmp1);

//	x-polarized beam
	double para=sinalpha*cosalpha*exp(fabs(DblA)-sinalpha*sinalpha*(1.0/sx/sx+1.0/sy/sy)/8.0);
	double paraA1=(NPAITAO[0]+cosalpha*NPAITAO[1])*para;
	double paraA2=(NPAITAO[0]-cosalpha*NPAITAO[1])*para;
	double paraB1=(NPAITAO[1]+cosalpha*NPAITAO[0])*para;
	double paraB2=(NPAITAO[1]-cosalpha*NPAITAO[0])*para;
	tmpA[0]=paraA1*Xmn1[0]+paraA2*Xmp1[0],tmpA[1]=paraA1*Xmn1[1]+paraA2*Xmp1[1];
	tmpB[0]=paraB1*Xmn1[0]+paraB2*Xmp1[0],tmpB[1]=paraB1*Xmn1[1]+paraB2*Xmp1[1];
/*
//	axisymetric beam
	double para=sinalpha*cosalpha*(1.0+cosalpha)*0.5*exp(fabs(DblA)-sinalpha*sinalpha*(1.0/sx/sx+1.0/sy/sy)/8.0);
	double paraA1=(NPAITAO[0]+NPAITAO[1])*para;
	double paraA2=(NPAITAO[0]-NPAITAO[1])*para;
	tmpA[0]=paraA1*Xmn1[0]+paraA2*Xmp1[0],tmpA[1]=paraA1*Xmn1[1]+paraA2*Xmp1[1];
	tmpB[0]=paraA1*Xmn1[0]-paraA2*Xmp1[0],tmpB[1]=paraA1*Xmn1[1]-paraA2*Xmp1[1];
*/
	para=Z0*cosalpha,tmp0[0]=cos(para),tmp0[1]=-sin(para);
	CplxMul(tmp0[0],tmp0[1],tmpA[0],tmpA[1],tmpA);	res[0]=tmpA[0],res[1]=tmpA[1];
	CplxMul(tmp0[0],tmp0[1],tmpB[0],tmpB[1],tmpB);	res[2]=tmpB[0],res[3]=tmpB[1];
	res[4]=Funcabs(res[0],res[1]),res[5]=Funcabs(res[2],res[3]);
	return;
}



//	**************************************************************************************************
//	**************************************************************************************************
void FuncXmBesselJB1(int m,double DblA,double Beta,double Ksi,double Xmn1[],double Xmp1[])
{
	int m1=m-1,m2=m+1;
//	**************************************************************************************************
	if(DblA==0.0)
	{
		double Jm1=Real_BesselJ_Downward_Single(m1,Beta);
		double Jm2=Real_BesselJ_Downward_Single(m2,Beta);
		Xmn1[0]=cos(m1*Ksi)*Jm1,Xmn1[1]=-sin(m1*Ksi)*Jm1;
		Xmp1[0]=cos(m2*Ksi)*Jm2,Xmp1[1]=-sin(m2*Ksi)*Jm2;
	}
//	**************************************************************************************************
	else
	{
		double tmp[2],absDblA=fabs(DblA),absBeta=fabs(Beta);
		int p,s,absp,abss;
		int pmax0=int(absBeta+10.55*pow(absBeta,0.327))+7;
		int smax0=int(8.312*sqrt(absDblA))+5;
		int s1=-(pmax0+m)/2,s2=(pmax0-m)/2;
		if(s1+smax0>0) s1=-smax0;
		if(s2<smax0) s2=smax0;
		int smax=((s1+s2)<0)?-s1:s2;
		int p1=abs(m+2*s1),p2=abs(m+2*s2);
		int pmax=(p1>p2)?p1:p2;

		double *EqBslI2A=(double*)malloc((smax+5)*sizeof(double));
		double *BslJBeta=(double*)malloc((pmax+5)*sizeof(double));
		Real_EqBesselI(smax+3,DblA,EqBslI2A);
		Real_BesselJ_Downward(pmax+3,Beta,BslJBeta);

		Xmn1[0]=Xmn1[1]=Xmp1[0]=Xmp1[1]=0.0;
		double pKsi,Vtmp;
		for(s=s1;s<=s2;s++)
		{
			abss=abs(s);
			p=m1+2*s,absp=abs(p),pKsi=p*Ksi;
			Vtmp=BslJBeta[absp]*EqBslI2A[abss];
			if(p<0 && absp%2==1) tmp[0]=-cos(pKsi)*Vtmp,tmp[1]=sin(pKsi)*Vtmp;
			else tmp[0]=cos(pKsi)*Vtmp,tmp[1]=-sin(pKsi)*Vtmp;
			if(abss%2==0) Xmn1[0]+=tmp[0],Xmn1[1]+=tmp[1];
			else Xmn1[0]-=tmp[0],Xmn1[1]-=tmp[1];

			p=m2+2*s,absp=abs(p),pKsi=p*Ksi;
			Vtmp=BslJBeta[absp]*EqBslI2A[abss];
			if(p<0 && absp%2==1) tmp[0]=-cos(pKsi)*Vtmp,tmp[1]=sin(pKsi)*Vtmp;
			else tmp[0]=cos(pKsi)*Vtmp,tmp[1]=-sin(pKsi)*Vtmp;
			if(abss%2==0) Xmp1[0]+=tmp[0],Xmp1[1]+=tmp[1];
			else Xmp1[0]-=tmp[0],Xmp1[1]-=tmp[1];
		}
		Xmn1[2]=Funcabs(Xmn1[0],Xmn1[1]),Xmp1[2]=Funcabs(Xmp1[0],Xmp1[1]);
		free(EqBslI2A),free(BslJBeta); EqBslI2A=BslJBeta=NULL;
	}
	return;
}


//	**************************************************************************************************
//	**************************************************************************************************
void FuncXmBesselJB(int m,double DblA,double Beta,double Ksi,double Xmn1[],double Xmp1[])
{
	int m1=m-1,m2=m+1;
//	**************************************************************************************************
	if(DblA==0.0)
	{
		double Jm1=Real_BesselJ_Downward_Single(m1,Beta);
		double Jm2=Real_BesselJ_Downward_Single(m2,Beta);
		Xmn1[0]=cos(m1*Ksi)*Jm1,Xmn1[1]=-sin(m1*Ksi)*Jm1;
		Xmp1[0]=cos(m2*Ksi)*Jm2,Xmp1[1]=-sin(m2*Ksi)*Jm2;
	}
//	**************************************************************************************************
	else
	{
		double tmp[2],absDblA=fabs(DblA),absBeta=fabs(Beta);
		int p,s,absp,abss,smax,pmax,s1,s2,p1,p2;
		if(absDblA<absBeta) pmax=int(absBeta+10.5*pow(absBeta,0.341)+5.5);
		else pmax=int(absBeta+10.5*pow(absBeta,0.341)+10.5);
		s1=-(pmax+m)/2,s2=(pmax-m)/2;
		if(m>0 && s2<20) s2=20;
		if(m<0 && s1>-20) s1=-20;
		smax=((s1+s2)<0)?-s1:s2;
		if((p1=m+2*s1)<0) p1=-p1;
		if((p2=m+2*s2)<0) p2=-p2;
		pmax=(p1>p2)?p1:p2;

		double *EqBslI2A=(double*)malloc((smax+5)*sizeof(double));
		double *BslJBeta=(double*)malloc((pmax+5)*sizeof(double));
		Real_EqBesselI(smax+3,DblA,EqBslI2A);
		Real_BesselJ_Downward(pmax+3,Beta,BslJBeta);

		Xmn1[0]=Xmn1[1]=Xmp1[0]=Xmp1[1]=0.0;
		double pKsi,Vtmp;
		for(s=s1;s<=s2;s++)
		{
			abss=abs(s);
			p=m1+2*s,absp=abs(p),pKsi=p*Ksi;
			Vtmp=BslJBeta[absp]*EqBslI2A[abss];
			if(p<0 && absp%2==1) tmp[0]=-cos(pKsi)*Vtmp,tmp[1]=sin(pKsi)*Vtmp;
			else tmp[0]=cos(pKsi)*Vtmp,tmp[1]=-sin(pKsi)*Vtmp;
			if(abss%2==0) Xmn1[0]+=tmp[0],Xmn1[1]+=tmp[1];
			else Xmn1[0]-=tmp[0],Xmn1[1]-=tmp[1];

			p=m2+2*s,absp=abs(p),pKsi=p*Ksi;
			Vtmp=BslJBeta[absp]*EqBslI2A[abss];
			if(p<0 && absp%2==1) tmp[0]=-cos(pKsi)*Vtmp,tmp[1]=sin(pKsi)*Vtmp;
			else tmp[0]=cos(pKsi)*Vtmp,tmp[1]=-sin(pKsi)*Vtmp;
			if(abss%2==0) Xmp1[0]+=tmp[0],Xmp1[1]+=tmp[1];
			else Xmp1[0]-=tmp[0],Xmp1[1]-=tmp[1];
		}
		Xmn1[2]=Funcabs(Xmn1[0],Xmn1[1]),Xmp1[2]=Funcabs(Xmp1[0],Xmp1[1]);
		free(EqBslI2A),free(BslJBeta); EqBslI2A=BslJBeta=NULL;
	}
	return;
}




//	**************************************************************************************************
//	**************************************************************************************************
int FuncRelaModifiedBessel(int m,double DblA[],double Beta[],double cosKsi[],double sinKsi[],double FIm[][2])
{
	int mn;
	double absDblA=Funcabs(DblA[0],DblA[1]);
	double absBeta=Funcabs(Beta[0],Beta[1]);
	double tmp0[2],tmp1[2],tmp2[2],tmp3[2];
	double PKsi[2]={cosKsi[0]+sinKsi[1],cosKsi[1]-sinKsi[0]};
	double NKsi[2]={cosKsi[0]-sinKsi[1],cosKsi[1]+sinKsi[0]};

	if(absDblA==0.0)
	{
		double resmn1[2],resm0[2],resmp1[2];
		CMPLE_EqBesselI_3Order(m,Beta[0],Beta[1],resmn1,resm0,resmp1);

		mn=m;
		if(mn>=0) CplxPow(mn,PKsi[0],PKsi[1],tmp0);
		else CplxPow(-mn,NKsi[0],NKsi[1],tmp0);
		CplxMul(resm0[0],resm0[1],tmp0[0],tmp0[1],tmp1);

		mn=m-1;
		if(mn>=0) CplxPow(mn,PKsi[0],PKsi[1],tmp0);
		else CplxPow(-mn,NKsi[0],NKsi[1],tmp0);
		CplxMul(resmn1[0],resmn1[1],tmp0[0],tmp0[1],tmp2);

		mn=m+1;
		if(mn>=0) CplxPow(mn,PKsi[0],PKsi[1],tmp0);
		else CplxPow(-mn,NKsi[0],NKsi[1],tmp0);
		CplxMul(resmp1[0],resmp1[1],tmp0[0],tmp0[1],tmp3);

		FIm[0][0]=tmp1[0],FIm[0][1]=tmp1[1];						//	I(m)
		FIm[1][0]=tmp3[0]+tmp2[0],FIm[1][1]=tmp3[1]+tmp2[1];		//	I(m+1)+I(m-1)
		FIm[2][0]=tmp3[0]-tmp2[0],FIm[2][1]=tmp3[1]-tmp2[1];		//	I(m+1)-I(m-1)
		return 1;
	}
	else
	{
		double V0=0.0,V12=0.0,Vtmp;
		int p,s,abss,abssp1,absp,smax,pmax,s1,s2,p1,p2,absm=abs(m);
		pmax=int(absBeta+10.5*pow(absBeta,0.341)+20.5);
		s1=-(pmax+m)/2,s2=(pmax-m)/2;
		if(m>0 && s2<10) s2=10;
		if(m<0 && s1>-10) s1=-10;
		smax=(s1+s2<0)?-s1:s2;
		if((p1=m+2*s1)<0) p1=-p1;
		if((p2=m+2*s2)<0) p2=-p2;
		pmax=(p1>p2)?p1:p2;

		double *reEqBslI2A=(double*)calloc((smax+5),sizeof(double));
		double *imEqBslI2A=(double*)calloc((smax+5),sizeof(double));
		double *reEqBslIBeta=(double*)calloc((pmax+5),sizeof(double));
		double *imEqBslIBeta=(double*)calloc((pmax+5),sizeof(double));
		CMPLE_EqBesselI(smax+3,DblA[0],DblA[1],reEqBslI2A,imEqBslI2A);
		CMPLE_EqBesselI(pmax+3,Beta[0],Beta[1],reEqBslIBeta,imEqBslIBeta);
//	=====================================================================================================
		FIm[0][0]=FIm[0][1]=FIm[1][0]=FIm[1][1]=FIm[2][0]=FIm[2][1]=0.0;
		for(s=s1;s<=s2;s++)
		{
			abss=abs(s);
			p=s+1,abssp1=abs(p);
//		I(m):
			p=m+2*s,absp=abs(p);
			if(p>=0) CplxPow(p,PKsi[0],PKsi[1],tmp0);
			else CplxPow(-p,NKsi[0],NKsi[1],tmp0);
			CplxMul(reEqBslIBeta[absp],imEqBslIBeta[absp],tmp0[0],tmp0[1],tmp0);
			CplxMul(reEqBslI2A[abss],imEqBslI2A[abss],tmp0[0],tmp0[1],tmp1);

//		I(m+1)+I(m-1) & I(m+1)-I(m-1):
			p=m+1+2*s,absp=abs(p);
			if(p>=0) CplxPow(p,PKsi[0],PKsi[1],tmp0);
			else CplxPow(-p,NKsi[0],NKsi[1],tmp0);
			CplxMul(reEqBslIBeta[absp],imEqBslIBeta[absp],tmp0[0],tmp0[1],tmp0);

			CplxMul(reEqBslI2A[abss]+reEqBslI2A[abssp1],imEqBslI2A[abss]+imEqBslI2A[abssp1],tmp0[0],tmp0[1],tmp2);
			CplxMul(reEqBslI2A[abss]-reEqBslI2A[abssp1],imEqBslI2A[abss]-imEqBslI2A[abssp1],tmp0[0],tmp0[1],tmp3);

			FIm[0][0]+=tmp1[0],FIm[0][1]+=tmp1[1];		//	I(m)
			FIm[1][0]+=tmp2[0],FIm[1][1]+=tmp2[1];		//	I(m+1)+I(m-1)
			FIm[2][0]+=tmp3[0],FIm[2][1]+=tmp3[1];		//	I(m+1)-I(m-1)

			if(V0<(Vtmp=Funcabs(tmp1[0],tmp1[1]))) V0=Vtmp;
			if(V12<(Vtmp=Funcabs(tmp2[0],tmp2[1]))) V12=Vtmp;
			if(V12<(Vtmp=Funcabs(tmp3[0],tmp3[1]))) V12=Vtmp;
		}
		free(reEqBslI2A),free(imEqBslI2A),free(reEqBslIBeta),free(imEqBslIBeta);
		reEqBslI2A=imEqBslI2A=reEqBslIBeta=imEqBslIBeta=NULL;

		if(Funcabs(FIm[0][0],FIm[0][1])<1.0e-14*V0)	 FIm[0][0]=FIm[0][1]=0.0;
		if(Funcabs(FIm[1][0],FIm[1][1])<1.0e-14*V12) FIm[1][0]=FIm[1][1]=0.0;
		if(Funcabs(FIm[2][0],FIm[2][1])<1.0e-14*V12) FIm[2][0]=FIm[2][1]=0.0;
		return 1;
	}
}



//	**************************************************************************************************
//	**************************************************************************************************
//	Normalized Associated Legendre Function: Pnm(costheta)
//	Ferrer definition
//	**************************************************************************************************
//	**************************************************************************************************
double Norm_AssLegFunc_Single(int n,int m, double theta)
{
	m=abs(m),n=abs(n);
	double miu=cos(theta),fn=n,fm=m;
	if(m==0)
	{
		if(n==0) return sqrt(0.5);
		else if(n==1) return (sqrt(1.5)*miu);
		else
		{
			double NPn0=sqrt(0.5),NPn1=sqrt(1.5)*miu,NPn2;
			for(int k=2;k<=n;k++)
			{
				double fk=k;
				NPn2=sqrt(4.0-1.0/fk/fk)*(miu*NPn1-NPn0/sqrt(4.0-1.0/(fk-1.0)/(fk-1.0)));
				NPn0=NPn1,NPn1=NPn2;
			}
			return NPn2;
		}
	}
	else
	{
		if(fabs(miu)==1.0 || n<m) return 0.0;
		else
		{
			int k;
			double Vmm=log(sin(theta))*m+0.5*log(fm+0.5);
			for(k=1;k<=m;k++) Vmm+=0.5*log(1.0-0.5/k);
			if(n==m) return (exp(Vmm));
			else
			{
				double Vtmp=Vmm,res0=0.0,res1=1.0,res2,Clg1e80=log(1e80);
				for(k=m+1;k<=n;k++)
				{
					double fk=k;
					res2=sqrt((4.0*fk*fk-1.0)/(fk*fk-fm*fm))*(miu*res1-sqrt(((fk-1.0)*(fk-1.0)-fm*fm)/(4.0*(fk-1.0)*(fk-1.0)-1.0))*res0);
					if(fabs(res2)>1e80) res0=res1*1e-80,res1=res2*1e-80,Vtmp+=Clg1e80;
					else if(fabs(res2)<1e-80) res0=res1*1e80,res1=res2*1e80,Vtmp-=Clg1e80;
					else res0=res1,res1=res2;
				}
				Vtmp=exp(Vtmp);
				return (res1*Vtmp);
			}
		}
	}
	return 0;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Normalized Associated Legendre Function: Pn0(costheta)
//	**************************************************************************************************
//	**************************************************************************************************
void Norm_AssLegFunc_n0(int nmax,double theta,double NALFn0[])
{
	double miu=cos(theta);
	NALFn0[0]=sqrt(0.5),NALFn0[1]=sqrt(1.5)*miu;
	for(int n=2;n<=nmax;n++) NALFn0[n]=(miu*NALFn0[n-1]-NALFn0[n-2]*(n-1.0)/sqrt((2.0*n-1.0)*(2.0*n-3.0)))*sqrt(4.0-1.0/n/n);
	return;
}



//	**************************************************************************************************
//	**************************************************************************************************
//	angular functions
//	**************************************************************************************************
//	**************************************************************************************************
void Norm_AngularFunc_Single(int n,int m,double theta,double NPAITAO[])
{
	double *NPAI=(double*)malloc((n+5)*sizeof(double));
	double *NTAO=(double*)malloc((n+5)*sizeof(double));
	Norm_AngularFunc_COL(n,theta,NPAI,NTAO);
	NPAITAO[0]=NPAI[m],NPAITAO[1]=NTAO[m];
	free(NPAI),free(NTAO),NPAI=NTAO=NULL;
	return;
}



//	**************************************************************************************************
//	**************************************************************************************************
//	fixed n,  m=0 -> n
//	**************************************************************************************************
//	**************************************************************************************************
void Norm_AngularFunc_COL(int n,double theta,double NPAI[],double NTAO[])
{
	int m,k;
	double sinx=sin(theta),cosx=cos(theta),fn=n;
	if(n==1)
	{
		NPAI[0]=0.0,NPAI[1]=0.5*sqrt(3.0);
		NTAO[0]=-sinx*sqrt(1.5),NTAO[1]=NPAI[1]*cosx;
	}
	else
	{
		if(sinx>0.0)
		{
			double cotx=cosx/sinx;
			NPAI[n+1]=NPAI[0]=0.0;
//	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			double Vnn=log(sinx)*(n-1)+log(n)+0.5*log(n+0.5);
			for(k=1;k<=n;k++) Vnn+=0.5*log(1.0-0.5/k);
			double Vtmp=Vnn,Clg1e80=log(1e80);
			NPAI[n]=1.0;
			for(m=n-1;m>=1;m--)
			{
				double fm=m;
				NPAI[m]=(2.0*cotx*NPAI[m+1]-NPAI[m+2]*sqrt((fn-fm-1)*(fn+fm+2))/(fm+2))*fm/sqrt((fn-fm)*(fn+fm+1));
				if(fabs(NPAI[m])>1e80)
				{
					for(k=m;k<=n;k++) NPAI[k]*=1e-80;
					Vtmp+=Clg1e80;
				}
				if(fabs(NPAI[m])<1e-80)
				{
					for(k=m;k<=n;k++) NPAI[k]*=1e80;
					Vtmp-=Clg1e80;
				}
			}
			if(Vtmp>709.0) printf("\n\t%le\n\n",Vtmp),getche();
			else Vtmp=exp(Vtmp);
			for(m=1;m<=n;m++) NPAI[m]*=Vtmp;
//	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			for(m=0;m<=n;m++) NTAO[m]=cosx*NPAI[m]-sinx*sqrt((n-m)*(n+m+1.0))*NPAI[m+1]/(m+1);
		}
		else
		{
			for(m=0;m<=n;m++) NPAI[m]=NTAO[m]=0.0;
			if(cosx>0.0 || (n+1)%2==0) NPAI[1]=sqrt(fn*(fn+0.5)*(fn+1.0))*0.5;
			else NPAI[1]=-sqrt(fn*(fn+0.5)*(fn+1.0))*0.5;
			NTAO[1]=(cosx>0.0)? NPAI[1]:-NPAI[1];
		}
	}
	return;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Bessel function: Jn(z)
//	**************************************************************************************************
//	**************************************************************************************************
void Real_BesselJ_ThreeOrder(int n, double z,double BslJ[])
{
	int absn=abs(n);
	double *BesselJ=(double*)malloc((absn+3)*sizeof(double));
	Real_BesselJ_Downward(absn+1,z,BesselJ);

	if(n==0) BslJ[0]=-BesselJ[1],BslJ[1]=BesselJ[0],BslJ[2]=BesselJ[1];
	else if(n>0) BslJ[0]=BesselJ[n-1],BslJ[1]=BesselJ[n],BslJ[2]=BesselJ[n+1];
	else
	{
		if(absn%2==0) BslJ[0]=-BesselJ[absn+1],BslJ[1]=BesselJ[absn],BslJ[2]=-BesselJ[absn-1];
		else BslJ[0]=BesselJ[absn+1],BslJ[1]=-BesselJ[absn],BslJ[2]=BesselJ[absn-1];
	}
	free(BesselJ); BesselJ=NULL;
	return;
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Jn(z): downward recurrence
//	**************************************************************************************************
//	**************************************************************************************************
void Real_BesselJ_Downward(int nmax, double z,double BesselJ[])
{
	int sign=0;
	if(nmax<0) nmax=-nmax,sign+=1;
	if(fabs(z)<=1e-30)
	{
		BesselJ[0]=1.0;
		for(int n=1;n<=nmax;n++) BesselJ[n]=0.0;
		return;
	}
	else
	{
//	***********************************************************************
		if(z<0.0) sign+=1,z=-z;
		int i,n,qplus;
		if(float(nmax)<z) qplus=int(z+10.5*pow(z,0.341)+2.5);
		else qplus=nmax;

		double BslJ2,BslJ1,BslJ0;
		do{
			qplus+=10,BslJ2=0.0,BslJ1=1.0;
			for(n=qplus;n>=nmax;n--)
			{
				BslJ0=2.0*(n+1)/z*BslJ1-BslJ2,BslJ2=BslJ1,BslJ1=BslJ0;
				if(BslJ0>1.0e18) break;
			}
		}
		while(BslJ0<1.0e18);

//	***********************************************************************

		double *BslJ=(double*)malloc((qplus+5)*sizeof(double));
		BslJ[qplus+2]=0.0,BslJ[qplus+1]=1e-10;
		for(n=qplus;n>=0;n--)
		{
			BslJ[n]=2.0*(n+1)/z*BslJ[n+1]-BslJ[n+2];
			double limit=fabs(BslJ[n]);
			if(limit>1e+80) {for(i=n;i<=qplus+1;i++) BslJ[i]/=limit;}
		}

//	***********************************************************************
		double tmp=0.0;
		for(n=1;n<=int(qplus/2);n++) tmp+=BslJ[2*n];
		tmp=tmp*2.0+BslJ[0];

		for(n=0;n<=nmax;n++) BesselJ[n]=((sign%2==1 && n%2==1)? (-BslJ[n]/tmp):(BslJ[n]/tmp));
		free(BslJ); BslJ=NULL;
		return;
	}
}


//	**************************************************************************************************
//	**************************************************************************************************
//	Jn(z) single order
//	**************************************************************************************************
//	**************************************************************************************************
double Real_BesselJ_Downward_Single(int n, double z)
{
	if(fabs(z)<1e-30) return ((n==0)?1.0:0.0);
	int sign=0;
	if(z<0.0) sign+=1,z=-z;
	if(n<0) sign+=1,n=-n;
//	***********************************************************************
	int i,k,qplus;
	if(float(n)<z) qplus=int(z+10.5*pow(z,0.341)+2.5);
	else qplus=n;

	double BslJ2,BslJ1,BslJ0;
	do {
		qplus+=10,BslJ2=0.0,BslJ1=1.0;
		for(k=qplus;k>=n;k--)
		{
			BslJ0=2.0*(k+1)/z*BslJ1-BslJ2,BslJ2=BslJ1,BslJ1=BslJ0;
			if(BslJ0>1.0e18) break;
		}
	}
	while(BslJ0<1.0e18);

//	***********************************************************************
	double *BslJ;
	BslJ=(double*)malloc((qplus+5)*sizeof(double));
	BslJ[qplus+2]=0.0,BslJ[qplus+1]=1e-10;
	for(k=qplus;k>=0;k--)
	{
		BslJ[k]=2.0*(k+1)/z*BslJ[k+1]-BslJ[k+2];
		double limit=fabs(BslJ[k]);
		if(limit>1e+80) {for(i=k;i<=qplus+1;i++) BslJ[i]/=limit;}
	}

//	***********************************************************************
	double tmp=0.0;
	for(k=1;k<=int(qplus/2);k++) tmp+=BslJ[2*k];
	tmp=tmp*2.0+BslJ[0];

	for(k=0;k<=n;k++) BslJ[k]/=tmp;
	double res=BslJ[n];
	free(BslJ); BslJ=NULL;
	if(sign%2==1 && n%2==1) res=-res;
	return res;
}



//	*****************************************************************************************
//	*****************************************************************************************
//	Modified Bessel function: 3 order
//	m-1��resmn1[]; 	m��resm0[]; 	m+1��resmp1[]
//	eqv��	In(z)*exp(-z), if rez>0
//			In(z)*exp(+z), if rez<0
//	*****************************************************************************************
//	*****************************************************************************************
void CMPLE_EqBesselI_3Order(int m,double rez,double imz,double resmn1[],double resm0[],double resmp1[])
{
//	***********************************************************************
	double absz=Funcabs(rez,imz);
	if(absz<1e-30)
	{
		resmn1[0]=resm0[0]=resmp1[0]=resmn1[1]=resm0[1]=resmp1[1]=0.0;
		if(m==0) resm0[0]=1.0;
		if(m==1) resmn1[0]=1.0;
		if(m==-1) resmp1[0]=1.0;
		return;
	}
//	***********************************************************************
	int n;
	if(m==0) n=1;
	else if(m>0) n=m+1;
	else n=1-m;
//	***********************************************************************
	double tmp[2],tmp1[2];
	int i,k,qplus;
	if(fabs(rez)>1.0 || n>int(absz+0.5)) qplus=n;
	else qplus=int(absz+7.5*pow(absz,0.34)+2.5);

	double BslI2[2],BslI1[2],BslI0[2],Vmax;
	do {
		qplus+=10,BslI2[0]=BslI2[1]=BslI1[1]=0.0,BslI1[0]=1.0;
		for(k=qplus;k>=n;k--)
		{
			CplxDiv(2.0*(k+1)*BslI1[0],2.0*(k+1)*BslI1[1],rez,imz,tmp);
			BslI0[0]=tmp[0]+BslI2[0],BslI0[1]=tmp[1]+BslI2[1];
			BslI2[0]=BslI1[0],BslI2[1]=BslI1[1];
			BslI1[0]=BslI0[0],BslI1[1]=BslI0[1];
			Vmax=Funcabs(BslI0[0],BslI0[1]);
			if(Vmax>1.0e15) break;
		}
	} while(Vmax<1.0e15);

	double *reBslI=(double*)malloc((qplus+5)*sizeof(double));
	double *imBslI=(double*)malloc((qplus+5)*sizeof(double));
//	***********************************************************************
	reBslI[qplus+2]=imBslI[qplus+2]=0.0;
	reBslI[qplus+1]=1e-10,imBslI[qplus+1]=0.0;
	for(k=qplus;k>=0;k--)
	{
		CplxDiv(2.0*(k+1)*reBslI[k+1],2.0*(k+1)*imBslI[k+1],rez,imz,tmp);
		reBslI[k]=tmp[0]+reBslI[k+2],imBslI[k]=tmp[1]+imBslI[k+2];
		double absBslI=Funcabs(reBslI[k],imBslI[k]);
		if(absBslI>1.0e80) {for(i=k;i<=qplus+1;i++) reBslI[i]/=absBslI,imBslI[i]/=absBslI;}
	}
//	***********************************************************************
	if(rez>=0.0)
	{
		tmp[0]=tmp[1]=0.0;
		for(k=1;k<=qplus;k++) tmp[0]+=reBslI[k],tmp[1]+=imBslI[k];
		tmp[0]=tmp[0]*2.0+reBslI[0],tmp[1]=tmp[1]*2.0+imBslI[0];
	}
	else
	{
		tmp[0]=tmp[1]=0.0;
		for(k=1;k<=qplus;k++)
		{
			if(k%2==0) tmp[0]+=reBslI[k],tmp[1]+=imBslI[k];
			else tmp[0]-=reBslI[k],tmp[1]-=imBslI[k];
		}
		tmp[0]=tmp[0]*2.0+reBslI[0],tmp[1]=tmp[1]*2.0+imBslI[0];
	}
	CplxDiv(1.0,0.0,tmp[0],tmp[1],tmp1);
	int kmin=(n>2)? (n-2):0;
	for(k=kmin;k<=n;k++)
	{
		CplxMul(reBslI[k],imBslI[k],tmp1[0],tmp1[1],tmp);
		reBslI[k]=tmp[0],imBslI[k]=tmp[1];
	}
	if(m==0) resmn1[0]=resmp1[0]=reBslI[1],resmn1[1]=resmp1[1]=imBslI[1],resm0[0]=reBslI[0],resm0[1]=imBslI[0];
	else if(m>0) resmn1[0]=reBslI[n-2],resmn1[1]=imBslI[n-2],resm0[0]=reBslI[n-1],resm0[1]=imBslI[n-1],resmp1[0]=reBslI[n],resmp1[1]=imBslI[n];
	else resmp1[0]=reBslI[n-2],resmp1[1]=imBslI[n-2],resm0[0]=reBslI[n-1],resm0[1]=imBslI[n-1],resmn1[0]=reBslI[n],resmn1[1]=imBslI[n];

	free(reBslI),free(imBslI); reBslI=imBslI=NULL;
	return;
}



//	*****************************************************************************************
//	*****************************************************************************************
//	Modified Bessel functions for complex variables
//	Bessel_I[n]=(reBesselI[n]+i*imBesselI[n])*exp(-z)	if rez>=0
//	Bessel_I[n]=(reBesselI[n]+i*imBesselI[n])*exp(z)	if rez<0
//	*****************************************************************************************
//	*****************************************************************************************
int CMPLE_EqBesselI(int nmax,double rez,double imz,double reBesselI[],double imBesselI[])
{
	if(nmax<0) nmax=-nmax;
	int n;
//	***********************************************************************
	double absz=Funcabs(rez,imz);
	if(absz<1e-30)
	{
		reBesselI[0]=1.0,imBesselI[0]=0.0;
		for(n=1;n<=nmax;n++) reBesselI[n]=imBesselI[n]=0.0;
	}
	else
	{
		double tmp1[2],tmp2[2];
		int k,qplus;
		if(fabs(rez)<1.0) qplus=int(absz+10.5*pow(absz,0.341)+2.5);
		if(qplus<nmax) qplus=nmax;

		double BslI2[2],BslI1[2],BslI0[2],Vmax;
		do {
			qplus+=10,BslI2[0]=BslI2[1]=BslI1[1]=0.0,BslI1[0]=1.0;
			for(n=qplus;n>=nmax;n--)
			{
				CplxDiv(2.0*(n+1)*BslI1[0],2.0*(n+1)*BslI1[1],rez,imz,tmp1);
				BslI0[0]=tmp1[0]+BslI2[0],BslI0[1]=tmp1[1]+BslI2[1];
				BslI2[0]=BslI1[0],BslI2[1]=BslI1[1];
				BslI1[0]=BslI0[0],BslI1[1]=BslI0[1];
				Vmax=Funcabs(BslI0[0],BslI0[1]);
				if(Vmax>1.0e15) break;
			}
		} while(Vmax<1.0e15);

		double *reBslI=(double*)malloc((qplus+5)*sizeof(double));
		double *imBslI=(double*)malloc((qplus+5)*sizeof(double));
//	***********************************************************************
		reBslI[qplus+2]=imBslI[qplus+2]=0.0;
		reBslI[qplus+1]=1e-10,imBslI[qplus+1]=0.0;
		for(n=qplus;n>=0;n--)
		{
			CplxDiv(2.0*(n+1)*reBslI[n+1],2.0*(n+1)*imBslI[n+1],rez,imz,tmp1);
			reBslI[n]=tmp1[0]+reBslI[n+2],imBslI[n]=tmp1[1]+imBslI[n+2];
			double absBslI=Funcabs(reBslI[n],imBslI[n]);
			if(absBslI>1.0e80) {for(k=n;k<=qplus+1;k++) reBslI[k]/=absBslI,imBslI[k]/=absBslI;}
		}
//	***********************************************************************
		if(rez>=0.0)
		{
			tmp1[0]=tmp1[1]=0.0;
			for(n=1;n<=qplus;n++) tmp1[0]+=reBslI[n],tmp1[1]+=imBslI[n];
			tmp1[0]=tmp1[0]*2.0+reBslI[0],tmp1[1]=tmp1[1]*2.0+imBslI[0];
		}
		else
		{
			tmp1[0]=tmp1[1]=0.0;
			for(k=1;k<=qplus/2;k++) n=2*k,tmp1[0]+=reBslI[n]-reBslI[n-1],tmp1[1]+=imBslI[n]-imBslI[n-1];
			tmp1[0]=tmp1[0]*2.0+reBslI[0],tmp1[1]=tmp1[1]*2.0+imBslI[0];
		}
		CplxDiv(1.0,0.0,tmp1[0],tmp1[1],tmp2);
		for(n=0;n<=nmax;n++)
		{
			CplxMul(reBslI[n],imBslI[n],tmp2[0],tmp2[1],tmp1);
			reBesselI[n]=tmp1[0],imBesselI[n]=tmp1[1];
		}
		free(reBslI),free(imBslI); reBslI=imBslI=NULL;
	}
	double V0=Funcabs(reBesselI[0],imBesselI[0])*1e-15;
	for(n=1;n<=nmax;n++) {if(Funcabs(reBesselI[n],imBesselI[n])<V0) {nmax=n;break;}}
	return nmax;
}


//	*****************************************************************************************
//	*****************************************************************************************
//	Modified Bessel functions for real variables
//	Bessel_I[n]=BesselI[n]*exp(-|z|)
//	*****************************************************************************************
//	*****************************************************************************************
int Real_EqBesselI(int nmax,double z,double BesselI[])
{
	int n;
	double absz=fabs(z);
//	***********************************************************************
	if(absz<1e-30)
	{
		BesselI[0]=1.0;
		for(n=1;n<=nmax;n++) BesselI[n]=0.0;
	}
	else
	{
//	***********************************************************************
		int k,qplus=nmax;
		double BslI2,BslI1,BslI0,tmp;
		do{
			qplus+=10,BslI2=0.0,BslI1=1.0;
			for(n=qplus;n>=nmax;n--)
			{
				BslI0=2.0*(n+1)/z*BslI1+BslI2;
				BslI2=BslI1,BslI1=BslI0;
				if(fabs(BslI0)>1.0e15) break;
			}
		} while(fabs(BslI0)<1.0e15);

		double *BslI=(double*)malloc((qplus+5)*sizeof(double));
//	***********************************************************************
		BslI[qplus+2]=0.0,BslI[qplus+1]=1e-10;
		for(n=qplus;n>=0;n--)
		{
			BslI[n]=2.0*(n+1)/z*BslI[n+1]+BslI[n+2];
			double tmp=fabs(BslI[n]);
			if(tmp>1.0e80) {for(k=n;k<=qplus+1;k++) BslI[k]/=tmp;}
		}
//	***********************************************************************
		if(z>=0.0) {for(tmp=0.0,n=1;n<=qplus;n++) tmp+=BslI[n];}
		else {for(tmp=0.0,k=1;k<=qplus/2;k++) n=2*k,tmp+=BslI[n]-BslI[n-1];}
		tmp=tmp*2.0+BslI[0];
		for(n=0;n<=nmax;n++) BesselI[n]=BslI[n]/tmp;
		free(BslI); BslI=NULL;
	}
	double V0=fabs(BesselI[0])*1e-15;
	for(n=1;n<=nmax;n++) if(fabs(BesselI[n])<V0) {nmax=n;break;}
	return nmax;
}


//	***************************************************************************************************
//	***************************************************************************************************
//	Riccati Bessel function of 1st kind
//	***************************************************************************************************
//	***************************************************************************************************
double Real_Riccati_BesselJ_Single(int n,double z)
{
	double absz=fabs(z),res;
	if(float(n)<absz)
	{
		double *RB1=(double*)malloc((n+5)*sizeof(double));
		RB1[0]=sin(z),RB1[1]=RB1[0]/z-cos(z);
		if(n>1) {for(int np=2;np<=n;np++) RB1[np]=(2*np-1.0)*RB1[np-1]/z-RB1[np-2];}
		res=RB1[n];
		free(RB1); RB1=NULL;
	}
	else
	{
		int np;
		int nmax=n+1;
		double *Dn1=(double*)malloc((nmax+5)*sizeof(double));
		double *RB1=(double*)malloc((nmax+5)*sizeof(double));
		double *tmp=(double*)malloc((nmax+5)*sizeof(double));
		LENTZ(nmax,z,0.0,tmp); Dn1[nmax]=tmp[0];
		for(np=nmax;np>=1;np--) tmp[np]=1.0/(Dn1[np]+np/z),Dn1[np-1]=np/z-tmp[np];
		RB1[0]=sin(z),RB1[1]=RB1[0]/z-cos(z);
		if(n>1) {for(np=2;np<=n;np++) RB1[np]=RB1[np-1]*tmp[np];}
		res=RB1[n];
		free(Dn1),free(RB1),free(tmp); Dn1=RB1=tmp=NULL;
	}
	return res;
}


//	***************************************************************************************************
//	***************************************************************************************************
//	Riccati Bessel function and its derivative, real variable
//	***************************************************************************************************
//	***************************************************************************************************
void Real_Riccati_BesselJ(int nmax,double z,double RB1[],double DRB1[])
//	***************************************************************************************************
//	***************************************************************************************************
{
	double absz=fabs(z);
	if(float(nmax)<absz)
	{
		RB1[0]=sin(z),DRB1[0]=cos(z);
		RB1[1]=RB1[0]/z-DRB1[0],DRB1[1]=RB1[0]-RB1[1]/z;
		if(nmax>1)
		{
			for(int np=2;np<=nmax;np++)
			{
				RB1[np]=(2*np-1.0)*RB1[np-1]/z-RB1[np-2];
				DRB1[np]=RB1[np-1]-np*RB1[np]/z;
			}
		}
	}
	else
	{
		double *Dn1=(double*)malloc((nmax+5)*sizeof(double));
		double *tmp=(double*)malloc((nmax+5)*sizeof(double));
		int n;
		LENTZ(nmax,z,0.0,tmp); Dn1[nmax]=tmp[0];
		for(n=nmax;n>=1;n--) tmp[n]=1.0/(Dn1[n]+n/z),Dn1[n-1]=n/z-tmp[n];
		RB1[0]=sin(z),RB1[1]=RB1[0]/z-cos(z),DRB1[1]=RB1[1]*Dn1[1];
		for(n=2;n<=nmax;n++) RB1[n]=RB1[n-1]*tmp[n],DRB1[n]=RB1[n]*Dn1[n];
		free(Dn1),free(tmp); Dn1=tmp=NULL;
	}
	return;
}




//	***************************************************************************************************
//	***************************************************************************************************
//	Riccati Bessel function and its derivative, complex variable
//	modified with exp(-imz)
//	***************************************************************************************************
//	***************************************************************************************************
void Complex_EqRiccati_BesselJ(int nmax,double rez,double imz,double reRB1[],double imRB1[],double reDRB1[],double imDRB1[])
{
	double absz=sqrt(rez*rez+imz*imz);
	int iabsz=int(absz+7.5*pow(absz,0.4))+20;
	int q=((nmax>iabsz)?nmax:iabsz);

	double *reDn1=(double*)malloc((q+5)*sizeof(double)),*imDn1=(double*)malloc((q+5)*sizeof(double));
	double *retmp=(double*)malloc((q+5)*sizeof(double)),*imtmp=(double*)malloc((q+5)*sizeof(double));

	double comx[2],comy[2];
	LENTZ(q,rez,imz,comx);
	reDn1[q]=comx[0],imDn1[q]=comx[1];

	if(imz==0.0)
	{
		int n;
		for(n=q;n>=1;n--) retmp[n]=1.0/(reDn1[n]+n/rez),reDn1[n-1]=n/rez-retmp[n];
		reRB1[0]=sin(rez),reRB1[1]=reRB1[0]/rez-cos(rez),reDRB1[1]=reRB1[1]*reDn1[1];
		imRB1[0]=imRB1[1]=imDRB1[1]=0.0;
		for(n=2;n<=nmax;n++)
		{
			reRB1[n]=reRB1[n-1]*retmp[n],reDRB1[n]=reRB1[n]*reDn1[n];
			imRB1[n]=imDRB1[n]=0.0;
		}
	}
	else
	{
		int n;
		for(n=q;n>=1;n--)
		{
			CplxDiv(double(n),0.0,rez,imz,comx);
			CplxDiv(1.0,0.0,comx[0]+reDn1[n],comx[1]+imDn1[n],comy);
			reDn1[n-1]=comx[0]-comy[0],imDn1[n-1]=comx[1]-comy[1];
			retmp[n]=comy[0],imtmp[n]=comy[1];
		}

		double expNdimz=exp(-2.0*imz);
		reRB1[0]=sin(rez)*0.5*(1.0+expNdimz);
		imRB1[0]=cos(rez)*0.5*(1.0-expNdimz);
		for(n=1;n<=nmax;n++)
		{
			CplxMul(reRB1[n-1],imRB1[n-1],retmp[n],imtmp[n],comx);
			reRB1[n]=comx[0],imRB1[n]=comx[1];
			CplxMul(reRB1[n],imRB1[n],reDn1[n],imDn1[n],comx);
			reDRB1[n]=comx[0],imDRB1[n]=comx[1];
		}
	}
	free(reDn1),free(imDn1),free(retmp),free(imtmp);
	reDn1=imDn1=retmp=imtmp=NULL;
	return;
}



//	***************************************************************************************************
//	***************************************************************************************************
//	Riccati Hankel function and its derivative, real variable
//	***************************************************************************************************
//	***************************************************************************************************
void Real_Ricatti_BesselH(int nmax,double z,double reRB3[],double imRB3[],double reDRB3[],double imDRB3[])
{
	double absz=fabs(z),sinz=sin(z),cosz=cos(z);
	reRB3[0]=sinz, reRB3[1]=sinz/z-cosz;
	reDRB3[0]=cosz,reDRB3[1]=sinz-reRB3[1]/z;
	imRB3[0]=-cosz,imRB3[1]=-(cosz/z+sinz);
	imDRB3[0]=sinz,imDRB3[1]=-cosz-imRB3[1]/z;
	if(nmax>1)
	{
		int n;
		if(float(nmax)<=absz)
			for(n=2;n<=nmax;n++) reRB3[n]=(2*n-1)*reRB3[n-1]/z-reRB3[n-2],reDRB3[n]=reRB3[n-1]-n*reRB3[n]/z;
		else
		{
			double *Dn1=(double*)malloc((nmax+5)*sizeof(double));
			LENTZ(nmax,z,0.0,Dn1); Dn1[nmax]=Dn1[0];
			for(n=nmax;n>=1;n--) Dn1[n-1]=n/z-1.0/(n/z+Dn1[n]);
			for(n=2;n<=nmax;n++) reRB3[n]=reRB3[n-1]/(n/z+Dn1[n]),reDRB3[n]=reRB3[n]*Dn1[n];
			free(Dn1); Dn1=NULL;
		}
		for(n=2;n<=nmax;n++) imRB3[n]=(2.0*n-1.0)*imRB3[n-1]/z-imRB3[n-2],imDRB3[n]=imRB3[n-1]-n*imRB3[n]/z;
	}
	return;
}


//	***************************************************************************************
//	***************************************************************************************
//	Lentz 
//	***************************************************************************************
//	***************************************************************************************
void LENTZ(int n,double rez,double imz,double res[])
{
	double abszs,reak,imak,reup,imup,relow,imlow,reFnk,imFnk,refact,imfact,tmp;
	abszs=rez*rez+imz*imz;
	reup=reFnk=(2.0*n+1.0)*rez/abszs,imup=imFnk=-(2.0*n+1.0)*imz/abszs;
	relow=imlow=0.0;

	for(int k=2;;k++)
	{
		tmp=(k%2==0)?(1.0-2.0*(n+k)):(2.0*(n+k)-1.0);
		reak=tmp*rez/abszs,imak=-tmp*imz/abszs;
		tmp=reup*reup+imup*imup,reup=reak+reup/tmp,imup=imak-imup/tmp;
		if(k==2) relow=reak,imlow=imak;
		else tmp=relow*relow+imlow*imlow,relow=reak+relow/tmp,imlow=imak-imlow/tmp;
		tmp=relow*relow+imlow*imlow;
		refact=(reup*relow+imup*imlow)/tmp,imfact=(imup*relow-reup*imlow)/tmp;
		res[0]=reFnk*refact-imFnk*imfact,res[1]=imFnk*refact+reFnk*imfact;
		if(fabs(refact*refact+imfact*imfact-1.0)<1.0e-30) break;
		else reFnk=res[0],imFnk=res[1];
	}
	res[0]-=rez/abszs*n,res[1]+=imz/abszs*n;
	return;
}



//	*****************************************************************************************
//	*****************************************************************************************
//	product of complex numbers (p+i*q)*(s+i*t)
//	*****************************************************************************************
//	*****************************************************************************************
void CplxMul(double p,double q,double s,double t,double x[])
{
	x[0]=p*s-q*t,x[1]=p*t+q*s;
	return;
}

//	*****************************************************************************************
//	*****************************************************************************************
//	(p+i*q)/(s+i*t)
//	*****************************************************************************************
//	*****************************************************************************************
void CplxDiv(double p,double q,double s,double t,double x[])
{
	double f1=(fabs(p)>fabs(q))? fabs(p):fabs(q);
	double f2=(fabs(s)>fabs(t))? fabs(s):fabs(t);
	if(f1>0.0 && f2>0.0)
	{
		p/=f1,q/=f1,s/=f2,t/=f2;
		double modl=s*s+t*t;
		x[0]=(p*s+q*t)/modl*(f1/f2),x[1]=(q*s-p*t)/modl*(f1/f2);
	}
	else x[0]=x[1]=0.0;
	return;
}

//	*****************************************************************************************
//	*****************************************************************************************
//	absolute of complex number
//	*****************************************************************************************
//	*****************************************************************************************
double Funcabs(double rez,double imz)
{
	if(rez<0.0) rez=-rez;
	if(imz<0.0) imz=-imz;
	if(rez>imz) return(rez*sqrt(1.0+pow(imz/rez,2)));
	else if(rez<imz) return(imz*sqrt(1.0+pow(rez/imz,2)));
	else return(rez*sqrt(2.0));
}

//	***************************************************************************************************
//	***************************************************************************************************
//	power of complex number
//	***************************************************************************************************
//	***************************************************************************************************
void CplxPow(int p,double reR,double imR,double res[])
{
	if(p==0) res[0]=1.0,res[1]=0.0;
	else if(reR==0.0 && imR==0.0) res[0]=res[1]=0.0;
	else
	{
		CplxLn(reR,imR,res);
		double a=res[0]*p,b=res[1]*p;
		res[0]=exp(a)*cos(b),res[1]=exp(a)*sin(b);
	}
	return;
}



//	***************************************************************************************************
//	***************************************************************************************************
//	log of complex number
//	***************************************************************************************************
//	***************************************************************************************************
void CplxLn(double a,double b,double res[])
{
	double mod=Funcabs(a,b);
	res[0]=log(mod);
	res[1]=(b>=0.0)? acos(a/mod):-acos(a/mod);
	return;
}
