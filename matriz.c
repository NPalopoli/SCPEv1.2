#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"


int matrix(int *Nsus_pos)

{
	
/*ORDEN DE AA ALFABETICO POR NOMBRE DE AA Ala, Arg,Asn,  ...*/

int orden[AA]={1,15,12,3,2,14,4,6,7,8,10,9,11,5,13,16,17,19,20,18};

FILE *outfile, *infile1;
FILE *outfile1, *outfile2;

	int i,j,a,b,c;
	double Nsus_sim[protsize][AA][AA];

	double Qij[AA][AA];
        double Qijref[AA][AA],freqref[AA]; 
        double Qijend[AA][AA],freqend[AA],Qijref_N[AA][AA];
        double Qij_rea[AA][AA];
        double Pij[AA][AA];
        double lambda_ref,lambda_rea,lambda,lambda_scpe;
        double alpha;

	double Sij[AA][AA],Sijref_N[AA][AA];
        double Sijend[AA][AA];
	double freq[AA],freq_rea[AA];
	double sum_total,sum[AA],suma, sum_q[AA];
	char filename[MAXNAME],filenamefreq[MAXNAME],filenameprob[MAXNAME];
        double lambda_chico=0.000001;

	double Norm(double freq[AA],double Qij[AA][AA]);
	double Sum_over_j(double N[AA]);
	int Print_Matrix(double N[AA][AA],FILE *out);
	
	infile1 = fopen ("qjones.dat", "r");
	sprintf(filenameprob,"%s.Prob",outmatrix_name);
	sprintf(filenamefreq,"%s.Freq",outmatrix_name);
     	outfile2 = fopen (filenameprob, "w");
        outfile1 = fopen (filenamefreq, "w");

printf("%s %s\n",filenameprob, filenamefreq);
        alpha=0.01;        
	

        /*Read reference Qij */
	for(i=0;i<AA;i++) { 
		for(j=0;j<AA;j++)  fscanf(infile1,"%lf ",&Qijref[i][j]);	
		fscanf(infile1,"\n");	
        }

        /*Read reference frequences*/
	for(i=0;i<AA;i++)  fscanf(infile1," %lf\n",&freqref[i]);

	
	/*Make Nsus symmetric*/ 
	for(c=0;c<protsize;c++)
		for(i=0;i<AA;i++)  
			for(j=0;j<AA;j++)  Nsus_sim[c][i][j]=0.5*(*(Nsus_pos + AA*AA*c +i*AA +j)+*(Nsus_pos + AA*AA*c +j*AA +i));

	for(c=0;c<protsize;c++) {
	/* Calculate Pij or Qij(i!=j)*/

		sum_total=0.;
		for(i=0;i<AA;i++) {
                         	sum[i]=0.;
				sum_q[i]=0.;
				sum[i]=Sum_over_j(Nsus_sim[c][i]);
				for(j=0;j<AA;j++) if(i!=j) {
						if(sum[i]!=0) Qij[i][j]=Nsus_sim[c][i][j]/sum[i];
					        	else Qij[i][j]=0.;
				} 
				sum_total+=sum[i];
				for(j=0;j<AA;j++) if(i!=j) sum_q[i]+=Qij[i][j]; 
				if(sum_q[i]!=0) Qij[i][i]=-sum_q[i];
                                       else Qij[i][i]=0.;
		}

	/*Calculate frecuencies*/
	for(i=0;i<AA;i++) if(sum_total!=0)freq[i]= sum[i]/sum_total;
					else freq[i]=0.;
	/*check SumFreq=1*/
				suma=0.;
				for(i=0;i<AA;i++) suma+=freq[i];
				if(suma!=1.) {
					for(i=0;i<AA;i++) if(suma!=0) freq[i]=freq[i]/suma;
				}

        /*Calcular lambda Qij ref*/
         lambda_ref=Norm(freqref,Qijref);

        /*Calcular lambda Qij scpe*/
	lambda_scpe=Norm(freq,Qij);


        /*NOrmalizar Qijref con respecto a Qij scpe*/
        for(i=0;i<AA;i++) for(j=0;j<AA;j++) 
                                      if(lambda_scpe!=0) Qijref_N[i][j]=(Qijref[i][j]*lambda_scpe)/lambda_ref;
                                      else Qijref_N[i][j]=Qijref[i][j]*lambda_chico; 

       /*Calcular Sij_ref*/
             for(i=0;i<AA;i++)
	         for(j=0;j<AA;j++) Sijref_N[i][j]= Qijref_N[i][j]/freqref[j];


       /*Aplicar seudocuentas a Freq[c]*/

	        for(i=0;i<AA;i++)
                    freqend[i]= (freq[i] + alpha*freqref[i])/ (1.0 + alpha);

         /*Calcular Sij */

             for(i=0;i<AA;i++)
	         for(j=0;j<AA;j++) Sij[i][j]= Qij[i][j]/freqend[j];


        /*Corregir Sij por seudocuentas*/
	for(i=0;i<AA;i++) for(j=0;j<AA;j++) 
               Sijend[i][j]= (Sij[i][j] + alpha*Sijref_N[i][j])/(alpha+1.0);


        /*Calculo Qijend a partir de Sijend*pj */
	for(i=0;i<AA;i++) for(j=0;j<AA;j++)
         Qijend[i][j]= Sijend[i][j] * freqend[j];

        /*renormalizo Qijend a la velocidad de Qij (scpe)*/
                               lambda=0.;
                               lambda=Norm(freqend,Qijend);
	                       for(i=0;i<AA;i++) for(j=0;j<AA;j++) 
                                      if(lambda!=0 || lambda_scpe!=0) Qijend[i][j]=(Qijend[i][j]*lambda_scpe)/lambda;
                                          else Qijend[i][j]=Qijref[i][j]*lambda_chico;

        /*Calculo Qijend a partir de Sijend*pj */
         for(i=0;i<AA;i++) for(j=0;j<AA;j++)
         Sijend[i][j]= Qijend[i][j]/freqend[j];

         /*OUTPUT SECTION*/
					

        /*Reordeno Qij Sij y Freq ya que HYPHY toma por orden alfabetico*/
        for(i=0;i<AA;i++) for(j=0;j<AA;j++) Sij[orden[i]-1][orden[j]-1]=Sijend[i][j];

        for(i=0;i<AA;i++) for(j=0;j<AA;j++) Qij_rea[orden[i]-1][orden[j]-1]=Qijend[i][j];

        for(i=0;i<AA;i++)     freq_rea[orden[i]-1]=freqend[i]; 



	/*Make HYPHY input*/
			sprintf(filename,"%s-%d.dat",outmatrix_name,c);
			outfile=fopen(filename,"w");
                        fprintf(outfile,"/*This file contains the substitution matrix for a given position or structural class obtained with scpe*/ \n\n");
                        fprintf(outfile,"  modelType = 0;\n");
                        fprintf(outfile,"  NICETY_LEVEL = 2;\n\n\n");
                        fprintf(outfile,"  %sp%d = {20,20};\n",Options.modelname,c+1);
			for(i=0;i<AA;i++) 
	             		for(j=i;j<AA;j++) {
                                 if(i!=j){
                                    fprintf(outfile," %sp%d[%d][%d] := a*%10.15f;\n",Options.modelname, c+1,i,j,Sij[i][j]);
                                    fprintf(outfile," %sp%d[%d][%d] := a*%10.15f;\n",Options.modelname, c+1,j,i,Sij[j][i]);
                                 }
			        }
			fprintf(outfile,"\n\n");
			fprintf(outfile,"freqs%s={20,1};\n",Options.modelname);
			for(i=0;i<AA;i++)  fprintf(outfile,"freqs%s[%d][0]=%10.8f; \n",Options.modelname,i,freq_rea[i]);
			fprintf(outfile,"\n\n");
			fprintf(outfile,"scpe%sp%d=0;\n",Options.modelname,c+1);
			fprintf(outfile,"Model %s = (%sp%d, freqs%s,1);\n",Options.modelname,Options.modelname,c+1,Options.modelname);
			fprintf(outfile,"FREQUENCY_SENSITIVE = 0;\n");
			fclose(outfile);



        /*Print Frecuencies*/

		for(i=0;i<AA;i++) fprintf(outfile1,"%f ",freq_rea[i]);
		fprintf(outfile1,"\n");
                   

        /*Print Probabilities*/

		for(i=0;i<AA;i++)sum[i]=0; 
		for(i=0;i<AA;i++) {
                                    for(j=0;j<AA;j++) {
                                              if(i!=j) sum[i]+= Qij_rea[i][j];
                                              Pij[i][j]= Qij_rea[i][j];
                                    }
                                    Pij[i][i]= 1. - sum[i];
               } 
              Print_Matrix(Pij,outfile2);



       }


     fclose(outfile1);
     fclose(outfile2);

   return(0);
}


/*funciones*/
						
double Sum_over_j(double N[AA])
{
						
		int i;
		double sum;
			
		sum=0.;
		for(i=0;i<AA;i++) sum+=N[i];
				
		return(sum);
							
}
						
int Print_Matrix(double N[AA][AA],FILE *out)
{
							
				
                int i,j;
						
		for(i=0;i<AA;i++) {
			for(j=0;j<AA;j++) fprintf(out,"%15.15f ",N[i][j]); 
			fprintf(out,"\n"); 
		}
					
		return(0);
}

double Norm(double freq[AA],double Qij[AA][AA])
{

		int i,j;
		double lambda;

                 
                lambda=0.;
		for(i=0;i<AA;i++) lambda+=Qij[i][i]*freq[i]; 

 

                return (lambda==0. ? 0. : -lambda);
}
/*
int Calculate_freq(double *pmatrix,double *freq)
{
  int i,j;
 
  for(i=0;i<AA;++i) *(freq+i)=(*(pmatrix+(i*AA+i)));
  return(0);
}*/

int Calculate_omega(double Nsus[AA][AA])

{

int i,j;
float sumN, sumS, omega,syn;


      sumN=sumS=syn=0.;

      for(i=0;i<AA;i++)
            for(j=0;j<AA;j++) if(i!=j) sumN += Nsus[i][j];
                                 else  sumS += Nsus[i][j];

      syn=(sumS - (sumN*(protsize-1)))/protsize;

      printf(" sumN %f sumS %f syn %f  omega %f\n",sumN,sumS,syn,sumN/(syn+sumN));


     return(0);
}
