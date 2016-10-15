/*UPO-SEVILLA iGEM 2016 TEAM. Information on how to use this software at: http://2016.igem.org/Team:UPO-Sevilla/Model  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main()
{	
	printf("\n");
	/*----------------------------------------------------------------------*/
	/* Declaration of all the parameters */
	double NF; /* Planktonic cell concentration in the reactor (mg/l) */
	double difNF; /* Differential of NF (mg/l) */
	double kmaxfree; /* Maximum planktonic cell growth rate */
	double NFVE; /* Input stream planktonic cell concentration */
	double KSF; /* Planktonic cell Monod constant */
	double KIfree; /* Planktonic cell inhibitory constant */
	double QE; /* Reactor input stream flow rate */
	double QS; /* Reactor output stream flow rate */
	double KF; /* Planktonic cell growth rate */
	double KB; /* Biofilm cell growth rate */
	double KA; /* Cell attachment rate */
	double KAbasalmin; /* Minimum basal cell attachment rate */
	double KAbasalmax; /* Maximum basal cell attachment rate */
	double KAbasal; /* Basal cell attachment rate */
	double KAinduced; /* Induced cell attachment rate */
	double GthA; /* Thresold glycerol concentration under which KAbasal equals KAbasalmin and over which KAbasal KAbasal equals KAbasalmax*/
	double KD; /* Cell detachment rate (1/hours)*/
	double KDbasalmin; /* Minimum basal cell detachment rate (1/hours)*/
	double KDbasalmax; /* Maximum basal cell detachment rate */
	double KDbasal; /* Basal cell detachment rate */
	double KDinduced; /* Induced cell detachment rate */
	double GthD; /* Thresold glycerol concentration under which 	KDbasal equals KDbasalmax and over which KDbasal equals KDbasalmin */
	double NB; /* Biofilm cell concentration in the reactor (mg/l) */
	double difNB; /* Differential of NB (mg/l))*/
	double NBVE; /* Input stream biofilm cell concentration (mg/l) */
	double NBVS; /* Output stream biofilm cell concentration (mg/l) */
	double kmaxbiofilm; /* Maximum biofilm cell growth rate */
	double KSB; /* Biofilm cell Monod constant (mg/l) */
	double KIbiofilm; /* Biofilm cell inhibitory constant (mg/l)*/
	double GVE; /* Input stream limiting substrate concentration */
	double G; /* Limiting substrate concentration (mg/l)*/
	double difG; /* Differential of G (mg)*/
	double t; /* Reaction time (hours)*/
	double dift;  /* Differential of t (hours)*/
	double tmax; /* Maximum reaction time (hours)*/
	double YF; /* Instantaneous ratio (difNF/dift) / difGNF , where difGNF represents the differential of limiting substrate consumed by the planktonic cells (strictly YF is a constant only in the stationary 	phase) */  
	double YB; /* Instantaneous ratio (difNB/dift) / difGNB , where difGNB represents the differential of limiting substrate consumed by the planktonic cells (strictly YB is a constant only in the 	stationary phase) */  
        double NFplusNB;
	FILE *data; /* For the file containing the data */
	FILE *results; /* For the file comprising the results */

	/*----------------------------------------------------------------------*/

	/* Data collection from a file, both of the initial values of the variable parameters and of the values of the constant parameters 	*/
	data=fopen("data_values.txt","rt");
	if(data==NULL)
	{
		printf("\ndata_values.txt error.");
	}
	while(fscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&NF,&difNF,&kmaxfree,&NFVE,&KSF,&KIfree,&QE,&QS,&KF,&KB,&KA,&KAbasalmin,&KAbasalmax,&KAbasal,&KAinduced,&GthA,&KD,&KDbasalmin,&KDbasalmax,&KDbasal,&KDinduced,&GthD,&NB,&difNB,&NBVE,&NBVS,&kmaxbiofilm,&KSB,&KIbiofilm,&GVE,&G,&difG,&t,&dift,&tmax,&YF,&YB)!=EOF)
	{
	}
	fclose(data);

	/*----------------------------------------------------------------------*/
           
        /* results_values.txt file preparation */
    	results=fopen("results_values.txt","at");
	if(results==NULL)
	{
		printf("\nresults_values.txt error");
	}


	/*----------------------------------------------------------------------*/

	/* Numerical integration to solve the differential equations.  */
	while (t<tmax)
	{	
		/* The first step is to figure out the value of variables which are necessary to know to calculate difNF, difNB and difG in the second step.*/
		KF=kmaxfree*G/(KSF+G);
		KB=kmaxbiofilm*G/(KSB+G);
		if(G>=GthA)
		{
			KAbasal=KAbasalmax; 
			/* Because GthA is the Thresold glycerol concentration underwhich KAbasal equals KAbasalmin and over which KAbasal KAbasal equals KAbasalmax */
		}
		if(G<GthA)
		{
			KAbasal=KAbasalmin;
			/* Because GthA is the Thresold glycerol concentration under which KAbasal equals KAbasalmin and over which KAbasal KAbasal equals KAbasalmax */
		}
		if(G<=GthD)
		{
			KDbasal=KDbasalmax;
			/* Because GthD is the Thresold glycerol concentration 	under which KDbasal equals KDbasalmax and over 	which KDbasal equals KDbasalmin */
		}
		if(G>GthD)
		{
			KDbasal=KDbasalmin;
			/* Because GthD is the Thresold glycerol concentration 			under 	which KDbasal equals KDbasalmax and over 				which KDbasal equals KDbasalmin */
		}
		KA=KAbasal+KAinduced;
		KD=KDbasal+KDinduced;
	
		/*The second step is to calculate difNF, difNB and difG. */
		difNF=(NFVE*QE-NF*QS+KF*NF-KA*NF+KD*NB)*dift;
		difNB=(NBVE*QE-NBVS*QS+KB*NB-KD*NB+KA*NF)*dift;
		if ((difNF>0)&&(difNB>0)) /* If there is both planktonic cell population growth and biofilm cell population growth */
		{
			difG=(GVE*QE-G*QS)*dift-(KF*NF*dift/YF)-(KB*NB*dift/YB);
		}
		if ((difNF<=0)&&(difNB>0)) /* If there is only biofilm cell population growth */
		{
			difG=(GVE*QE-G*QS)*dift-(KB*NB*dift/YB);
		}
		if ((difNF>0)&&(difNB<=0)) /* If there is only planktonic cell population growth */
		{
			difG=(GVE*QE-G*QS)*dift-(KF*NF*dift/YF);
		}
		if ((difNF<=0)&&(difNB<=0)) /* If there is neither planktonic cell population growth nor biofilm cell population growth */
		{
			difG=(GVE*QE-G*QS)*dift;
		}

		/*The third step consists in calculating the new value of NF, NB, G and t*/
		NF=NF+difNF;
		NB=NB+difNB;
		G=G+difG;
		NFplusNB=NF+NB;
		if(NF<0) /* Which would have no physical sense */
		{
			NF=0;
		}
		if(NB<0) /* Which would have no physical sense */
		{
			NB=0;
		}
		if(G<0) /* Which would have no physical sense */
		{
			G=0;
		}
		t=t+dift;

		/* In the fourth step the results are printed in the results_value.txt file*/
		fprintf(results,"%lf %lf %lf \n",t,NF,NB);
		
	}

	/* END */
	
	return 0;