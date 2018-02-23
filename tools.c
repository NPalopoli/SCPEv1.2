#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"

int ReadAlignment(FILE *ali_input)

{

int i,j,num;
char line[MAXSTR],str[MAXSTR];
int posit;
                i=0;
                while (fgets(line, MAXSTR, ali_input) != NULL ) {
                   if(line[0]=='>'){
                                 sscanf(line,">%s\n",ali[i].name);
                                 while (posit=ftell(ali_input),fgets(line, MAXSTR, ali_input) != NULL) {
                                    if(line[0]!='>' && isprint(line[0])){
                                               sscanf(line,"%s\n",str);
                                               strcat(ali[i].seq,str);
                                     }
                                    else break;

                                  }
                                  fseek(ali_input,posit,SEEK_SET);
                                  i++;

                   }
                }

        num=i;
       /*first sequence in alignment should be reference sequence*/
       for(i=0;i<num;i++)
            for(j=0;j<strlen(ali[0].seq);j++)
                        if(ali[i].seq[j]=='-') ali[i].seq[j]=ali[0].seq[j];

       for(j=0;j<num;j++) printf("%s\n%s\n",ali[j].name,ali[j].seq);

        return(i);

}

int Readseq(FILE *input,char *seq)

{
  int count;
  char name[MAXNAME],str[MAXSTR],line[MAXSTR];
				 
  count=0;
  fscanf(input,">%[^\n]\n",name);
  while (fgets(line, MAXSTR, input) != NULL /*, line[0]!='\n'*/ ) {
    sscanf(line,"%s\n",str);
    strcat(seq,str);
  }

  printf("Input Secuence\n%s\n%s\n",name,seq);
				 
  return(0);
}

/* Read PDB file */

int Readpdb(FILE *pdbin)

{

  FILE *contactfile;

  typedef struct {

    char atomgroup[4];
    char resname[4];
    char chain[2];
    int resnum;
    float coord[3];
    char flag;
    int sidechain;
    float vw;
    int sequential;
  }Atom;

  int i,j, natoms;
  int protlen,thisresidue;
  int count,resnum,*p;
  char line[MAXSTR], field[9];
  Atom atom[MAXATOMS] ;
  char title[MAXNAME];
  char *pline;
  char pdbseq[MAXPROTLEN];
  char pdbseqtotal[MAXPROTLEN];
  float xsum,ysum,zsum;
  float dist;
  int sum,res,flag,sumQ;
  int cont;

  
  p=&resnum;
  

  i=0;
  while(fgets(line,MAXSTR,pdbin)!=NULL){
    if(strstr(line,"ENDMDL")) break;
    if(strncmp(line,"TITLE",5)==0 && strlen(title)==0){
      sscanf(line,"TITLE %[^\n]",title);
      fprintf(logfile,"PDB TITLE %s\n",title);
    }

    if(strncmp(line,"ATOM",4)==0){

      strncpy(atom[i].atomgroup,line+13,3);
      atom[i].atomgroup[3]='\0';
      strncpy(atom[i].resname,line+17,3);
      atom[i].resname[3]='\0';

      strncpy(atom[i].chain,line+21,1);
      atom[i].chain[1]='\0';
      if(isspace(atom[i].chain[0]) || Options.chain=='_') {
                                          atom[i].chain[0]='u';
                                          Options.chain='u';
      }

     /* if(Options.chain=='a') {
                         atom[i].chain[0]='a';
                         Options.multimer=1;
      }*/
                     

      field[8]='\0';
      strncpy(field,line+23,6);
      atom[i].resnum=atoi(field);
										  
      /*Detect multiple comformations and choose one*/
      atom[i].flag=*(line+16);
      if(!strcmp(atom[i].atomgroup,atom[i-1].atomgroup)){
/*                  fprintf(stdout,"Residue %d has multiple comformations at %s atom. %c form is stored\n",atom[i].resnum,atom[i].atomgroup,atom[i-1].flag);*/
                  continue;
      }

      strncpy(field,line+30,8);
      atom[i].coord[0]=atof(field);
      strncpy(field,line+38,8);
      atom[i].coord[1]=atof(field);
      strncpy(field,line+46,8);
      atom[i].coord[2]=atof(field);

      /*Detect CA and sidechain atoms*/
      if(strncmp(atom[i].atomgroup,"N ",2)==0 || strncmp(atom[i].atomgroup,"O ",2)==0 || strncmp(atom[i].atomgroup,"C ",2)==0) atom[i].sidechain=0;
       if(atom[i].atomgroup[1]!=' ') atom[i].sidechain=1;


      /*Assign VdW radii to each atom*/
      if(atom[i].sidechain==1){
       if(strncmp(atom[i].atomgroup,"N",1)==0) atom[i].vw=1.70;
       if(strncmp(atom[i].atomgroup,"C",1)==0) atom[i].vw=1.80;
       if(strncmp(atom[i].atomgroup,"O",1)==0) atom[i].vw=1.40;
       if(strncmp(atom[i].atomgroup,"S",1)==0) atom[i].vw=2.00;
       if(strncmp(atom[i].atomgroup,"P",1)==0) atom[i].vw=2.00;
      }

					   
      i++;
      if(i==MAXATOMS-1){
	fprintf(stderr,"MAXATOMS limit has been reachead!!\n"); 
	exit(1);
      }
    }
  }

  natoms=i;


  /*extract and print selected protein sequence in PDB*/
  for(i=0,j=0;i<natoms;i++) 
    if(strstr(atom[i].atomgroup,"CA") && atom[i].chain[0]==Options.chain){
      Seq_to_number("three",atom[i].resname,p); 
      pdbseq[j]=aa_oneletter[*p];
      j++;
    }
    pdbseq[j]='\0';

   protsize=strlen(pdbseq);
   strcpy(seqprot,pdbseq);
   printf(">Protein chain %c with len %d\n%s\n",Options.chain,strlen(pdbseq),pdbseq);

   fprintf(logfile,">Protein chain %c with len %d\n%s\n",Options.chain,strlen(pdbseq),pdbseq);
  

/*extract ALL protein sequences in PDB*/
  for(i=0,j=0;i<natoms;i++) 
    if(strstr(atom[i].atomgroup,"CA")){
      Seq_to_number("three",atom[i].resname,p); 
      pdbseqtotal[j]=aa_oneletter[*p];
      j++;
    }
    pdbseqtotal[j]='\0';

   total_length=strlen(pdbseqtotal);
   strcpy(seqprot_total,pdbseqtotal);
   fprintf(logfile,"Total length of the PDB file is %d\n",total_length);
   fprintf(logfile,"%s\n",seqprot_total);


/*Assign sequential numbers to all residues in PDB. */

                   thisresidue=i=0;
                   count=1;
                   while(i<natoms){
                           thisresidue=atom[i].resnum;
                           while(thisresidue==atom[i].resnum){
                                   atom[i].sequential=count;
                                   i++;
                           }
                           count++;
                   }

 /*Calculate cmat*/

 if( (cmat=(int *)malloc((protsize*total_length)*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }

 if( (cmatQ=(int *)malloc((protsize*total_length)*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }

 for(i=0;i<protsize;++i) for(j=0;j<total_length;++j) *(cmat+(i*total_length+j))=0;
 for(i=0;i<protsize;++i) for(j=0;j<total_length;++j) *(cmatQ+(i*total_length+j))=0;

/*Extract  initial  number for selected chain*/

                   i=0;
                   while(i<natoms){
                           if(strstr(atom[i].atomgroup,"CA") && atom[i].chain[0]==Options.chain){
                                initial= atom[i].sequential;
                                break;   
                           }
                           i++;
                   }

printf("initial = %d \n",initial);

/*Monomer cmat*/
    i=thisresidue=0;
    while(i<natoms){
    thisresidue=atom[i].sequential;
    if(atom[i].sidechain==1 && atom[i].chain[0]== Options.chain){
       for(j=0;j<natoms;j++){
            if( atom[j].chain[0]==Options.chain){        
             if(atom[j].sequential < thisresidue-1 ||   atom[j].sequential > thisresidue+1) {
                   if(atom[j].sidechain==1){
                       dist=(atom[j].coord[0]-atom[i].coord[0])*(atom[j].coord[0]-atom[i].coord[0])+(atom[j].coord[1]-atom[i].coord[1])*(atom[j].coord[1]-atom[i].coord[1])+(atom[j].coord[2]-atom[i].coord[2])*(atom[j].coord[2]-atom[i].coord[2]);
                       if((sqrt(dist) - atom[i].vw - atom[j].vw) < 1.00){ 
/*                      printf("%d %d %c %d %d %s %s %f %f %f %f \n ",atom[i].sequential,atom[j].sequential,atom[j].chain[0],i+1,j+1,atom[i].atomgroup,atom[j].atomgroup, sqrt(dist),atom[i].vw,atom[j].vw,(sqrt(dist)-atom[i].vw-atom[j].vw));*/
	                                    *(cmat+((atom[i].sequential-initial)*total_length+atom[j].sequential-1))=1;

                       }
                   }
             }
          }
      }
    }
   i++;
   }

/*Quaternary cmat*/

    i=thisresidue=0;
    while(i<natoms){
    thisresidue=atom[i].sequential;
    if(atom[i].sidechain==1 && atom[i].chain[0]== Options.chain){
       for(j=0;j<natoms;j++){
           if(atom[j].chain[0]!=Options.chain){
             if(atom[j].sequential < thisresidue-1 || atom[j].sequential > thisresidue+1) {
                   if(atom[j].sidechain==1){
                       dist=(atom[j].coord[0]-atom[i].coord[0])*(atom[j].coord[0]-atom[i].coord[0])+(atom[j].coord[1]-atom[i].coord[1])*(atom[j].coord[1]-atom[i].coord[1])+(atom[j].coord[2]-atom[i].coord[2])*(atom[j].coord[2]-atom[i].coord[2]);
                       if((sqrt(dist) - atom[i].vw - atom[j].vw) < 1.00){ 
/*                     printf("%d %d %c %d %d %s %s %f %f %f %f \n",atom[i].sequential,atom[j].sequential,atom[j].chain[0],i+1,j+1,atom[i].atomgroup,atom[j].atomgroup, sqrt(dist),atom[i].vw,atom[j].vw,(sqrt(dist)-atom[i].vw-atom[j].vw));*/
                                      /* if(atom[j].chain[0]!=atom[i].chain[0]) printf(" contacto cuaternario en %d y %d\n",atom[i].sequential,res);*/
	                                    *(cmatQ+((atom[i].sequential-initial)*total_length+atom[j].sequential-1))=1;

                       }
                   }
                }
             }
       }
    }
    i++;
   }

  if((contactfile=fopen("Contact-information.dat","w"))==NULL) {
      fprintf(stderr,"Failed Output file for Quaternary contacts\n");
      exit(1);
    }

  for(i=0;i<protsize;++i)for(j=0;j<MAXCONTACTS;j++) cmatvector[i][j]=0;
  for(i=0;i<protsize;++i)for(j=0;j<MAXCONTACTS;j++) cmatQvector[i][j]=0;

  for(i=0;i<protsize;++i) {
                   cont=0; 
                   for(j=0;j<total_length;++j){
                          if(*(cmat+(i*total_length+j))==1){
                                 cmatvector[i][cont]=j;
                                 cont++;
                          }
                   }
                  if (cont>MAXCONTACTS){
                       fprintf(stderr,"increase MAXCONTACTS value\n");
                       exit(1);
                  }
  }

        
  for(i=0;i<protsize;++i) {
                   cont=0; 
                   for(j=0;j<total_length;++j){
                          if(*(cmatQ+(i*total_length+j))==1){
                                 cmatQvector[i][cont]=j;
                                 cont++;
                          }
                   }
                   if (cont>MAXCONTACTS){
                        fprintf(stderr,"increase MAXCONTACTS value\n");
                        exit(1);
                   }
 }


            fprintf(contactfile,"Tertiary contact map\n");
  for(i=0;i<protsize;++i) {
            sum=0;
            for(j=0;j<total_length;++j) sum+= *(cmat+(i*total_length+j));
            if(sum==0) fprintf(contactfile,"Warning!. Position %d has 0 contacts\n",i+1);
  }

            fprintf(contactfile,"Quaternary contact map\n");
  for(i=0;i<protsize;++i) {
            sum=0;
            for(j=0;j<total_length;++j) sum+= (*(cmatQ+(i*total_length+j))+ *(cmat+(i*total_length+j)));
            if(sum==0) fprintf(contactfile,"Warning!. Position %d has 0 contacts\n",i+1);
  }

            fprintf(contactfile,"Terciary \n");
  for(i=0;i<protsize;++i) {
            for(j=0;j<total_length;++j) fprintf(contactfile,"%d ",*(cmat+(i*total_length+j)));
            fprintf(contactfile," \n");
  }


            fprintf(contactfile,"Quaternary \n");
  for(i=0;i<protsize;++i) {
            for(j=0;j<total_length;++j) fprintf(contactfile,"%d ",(*(cmatQ+(i*total_length+j))+ *(cmat+(i*total_length+j))));
            fprintf(contactfile," \n");
  }


            fprintf(contactfile,"Tertiary Total contacts  Quaternary Total contacts    \n");

  for(i=0;i<protsize;++i) {
            sum=0;
            for(j=0;j<total_length;++j) sum+= *(cmat+(i*total_length+j));
            fprintf(contactfile,"%d ",sum);
            sumQ=0;
            for(j=0;j<total_length;++j) sumQ+= (*(cmatQ+(i*total_length+j))+ *(cmat+(i*total_length+j)));
            fprintf(contactfile,"          %d",sumQ);
            fprintf(contactfile,"  %d ", (sumQ - sum));
            fprintf(contactfile," \n");
  }
  return(0);

} 
				
int set_eijmat(float eijmat[])
{
  int iaa,jaa,kij;
  FILE *fp;
		  
  if((fp=fopen("ematrix.dat","r"))==NULL) {
      fprintf(stderr,"Contact potential matrix file  requiered\n");
      exit(1);
    }

      printf("\n using external ematrix file to read pair potential \n");
			  
      for(iaa=0;iaa<AA;iaa++){
	for(jaa=0;jaa<iaa+1;jaa++){
	  fscanf(fp,"%f",&eijmat[iaa*AA+jaa]);
	  eijmat[jaa*AA+iaa]=eijmat[iaa*AA+jaa];
	}
      }

  return(0);
		  
}	        
int Readq (FILE *qmutfile, double *qmut)

{
  int i,j;
  double sum;

  for(i=0;i<AA;++i){
    for(j=0;j<AA;++j){
      fscanf(qmutfile,"%lf ",qmut+(i*AA+j));
    }
    fgetc(qmutfile);
  }

  		 /*   printf("\nQmatrix inicial\n");
		    for(i=0;i<AA;++i){
		    for(j=0;j<AA;++j) printf("%15.15f ",*(qmut+(i*AA+j)));
		    printf("\n");
		    }*/
  return(0);
}

int Setp_cumulative(double *pmatrix)
{
  int i,j;
  double sum,*P,*Q;


  /*check P rows sum 1.00*/
  for(i=0;i<AA;++i){
    sum=0.;
    for(j=0;j<AA;++j) sum+=(*(pmatrix+(i*AA+j)));
    if(sum!=1.0) for(j=0;j<AA;++j) *(pmatrix+(i*AA+j))=((*(pmatrix+(i*AA+j)))/sum);
  }

  P=Q=pmatrix;
  for(i=0;i<AA;++i){
    for(j=1;j<AA;++j) {
      P++; 
      *P+=(*Q);
      Q++;
    }
    P++;Q++;
  } 

  	 /*   printf("P matrix cumulative\n");
	    for(i=0;i<AA;++i){
	    for(j=0;j<AA;++j) printf("%15.15f ",*(pmatrix+(i*AA+j)));
	    printf("\n");
	    }gutavo*/

  return(0);
}

int Setrate(double *qmatrix, double *freq)
{
  double lambda;
  int i,j;
       
  lambda=0.;
  lambda=Calculate_rate(qmatrix,freq);


  for(i=0;i<AA;++i) for(j=0;j<AA;++j) *(qmatrix+(i*AA+j))=(*(qmatrix+(i*AA+j)))/lambda;
  lambda=Calculate_rate(qmatrix,freq);

  return(0);
}
int Calculate_freq(double *pmatrix,double *freq)
{
  int i,j;
  double sum;

  for(i=0;i<AA;++i) *(freq+i)=(*(pmatrix+(i*AA+i)));
/*  sum=0.;
  for(j=0;j<AA;++j) {
          printf("freq %15.15f\n",*(freq+j));
          sum+=*(freq+j);
  }
          printf("sum =  %15.15f\n",sum);*/
  return(0);
}


double Calculate_rate(double *qmatrix,double *freq)

{
  int i,j;
  double lambda;
        
  lambda=0.;
  for(i=0;i<AA;++i) lambda+=(*(qmatrix+(i*AA+i)))*(*(freq+i));

  return(-lambda);
}



int Seq_to_number(char *code,char *seq, int *seq_num)

{
  if(!strcmp(code,"three")){
    while(*seq){
      if(strcmp("ALA",seq)==0) *seq_num=0;
      if(strcmp("ARG",seq)==0) *seq_num=1;
      if(strcmp("ASN",seq)==0) *seq_num=2;
      if(strcmp("ASP",seq)==0) *seq_num=3;
      if(strcmp("CYS",seq)==0) *seq_num=4;
      if(strcmp("GLN",seq)==0) *seq_num=5;
      if(strcmp("GLU",seq)==0) *seq_num=6;
      if(strcmp("GLY",seq)==0) *seq_num=7;
      if(strcmp("HIS",seq)==0) *seq_num=8;
      if(strcmp("ILE",seq)==0) *seq_num=9;
      if(strcmp("LEU",seq)==0) *seq_num=10;
      if(strcmp("LYS",seq)==0) *seq_num=11;
      if(strcmp("MET",seq)==0) *seq_num=12;
      if(strcmp("PHE",seq)==0) *seq_num=13;
      if(strcmp("PRO",seq)==0) *seq_num=14;
      if(strcmp("SER",seq)==0) *seq_num=15;
      if(strcmp("THR",seq)==0) *seq_num=16;
      if(strcmp("TRP",seq)==0) *seq_num=17;
      if(strcmp("TYR",seq)==0) *seq_num=18;
      if(strcmp("VAL",seq)==0) *seq_num=19;
      seq++;
      seq_num++;
    }
  }
  else
    while(*seq){
      if(islower(*seq)) *seq_num=20; /* position not structurally conserved */
      else
	switch(*seq){
	case 'A':
	  *seq_num=0;
	  break;
	case 'R':
	  *seq_num=1;
	  break;
	case 'N':
	  *seq_num=2;
	  break;
	case 'D':
	  *seq_num=3;
	  break;
	case 'C':
	  *seq_num=4;
	  break;
	case 'Q':
	  *seq_num=5;
	  break;
	case 'E':
	  *seq_num=6;
	  break;
	case 'G':
	  *seq_num=7;
	  break;
	case 'H':
	  *seq_num=8;
	  break;
	case 'I':
	  *seq_num=9;
	  break;
	case 'L':
	  *seq_num=10;
	  break;
	case 'K':
	  *seq_num=11;
	  break;
	case 'M':
	  *seq_num=12;
	  break;
	case 'F':
	  *seq_num=13;
	  break;
	case 'P':
	  *seq_num=14;
	  break;
	case 'S':
	  *seq_num=15;
	  break;
	case 'T':
	  *seq_num=16;
	  break;
	case 'W':
	  *seq_num=17;
	  break;
	case 'Y':
	  *seq_num=18;
	  break;
	case 'V':
	  *seq_num=19;
	  break;
	case '~':
	  *seq_num=20;
	  break;
	case '-':
	  *seq_num=20;
	  break;
	case '.':
	  *seq_num=20;
	  break;
	}
      seq++;
      seq_num++;
    }
				
				
  return(0);
}

int Number_to_seq(int *seqnum,char *seq)
{
  int cont;

  cont=0;
  while(cont<protsize){
    switch(*(seqnum+cont)){
    case 0:
      *seq='A';
      break;
    case 1:
      *seq='R';
      break;
    case 2:
      *seq='N';
      break;
    case 3:
      *seq='D';
      break;
    case 4:
      *seq='C';
      break;
    case 5:
      *seq='Q';
      break;
    case 6:
      *seq='E';
      break;
    case 7:
      *seq='G';
      break;
    case 8:
      *seq='H';
      break;
    case 9:
      *seq='I';
      break;
    case 10:
      *seq='L';
      break;
    case 11:
      *seq='K';
      break;
    case 12:
      *seq='M';
      break;
    case 13:
      *seq='F';
      break;
    case 14:
      *seq='P';
      break;
    case 15:
      *seq='S';
      break;
    case 16:
      *seq='T';
      break;
    case 17:
      *seq='W';
      break;
    case 18:
      *seq='Y';
      break;
    case 19:
      *seq='V';
      break;
    }
    cont++;
    seq++;
  }

  seq++;
  *seq='\0';
  return(0);
}



int MutateSeq(int *seq,double *pmatrix,int *nsyn)
{
  int i,aa,*P;
  int mut,p,cont;
  double r;

  P=seq;
  mut=0;
  *nsyn=0;

  if(!Options.screening){
   for (i=0; i<protsize; i++) {
    aa=*(P+i);
    *(P+i)=MakeMut(pmatrix+((*(P+i))*AA));
    if(aa!=(*(P+i))) {
      mut++;
       /*printf("mutation in pos=%d aa=%d *P=%d \n ",i, aa, *(P+i));*/
    }
   }
  }
 else{
       mut=0;
       r=rndu();
       if(r!=1.0) p=r*protsize;
         else p=0;
       cont=0;
       if(nsyn_persite[p]/Options.runs >= print_parameter){

             while(nsyn_persite[p]/Options.runs > print_parameter){

                      r=rndu();
                      if(r!=1.0) p=r*protsize;
                          else p=0;
                      vector[p]=1;

                    cont=0;
                    for(i=0;i<protsize;i++) {
                              if(vector[i]!=0){
                                     cont++;
                               }
                    }

                     if(cont==protsize) {
                                  average=1.;
                                  return(1);
                    }
              }
/*       for(i=0;i<protsize;i++)printf("%d ",vector[i]);
       printf("average %f cont %f pp=%f p=%d nsyn %f %f\n",average,number_of,print_parameter,p+1,nsyn_persite[p]/Options.runs,beta);*/

        }
/*       printf("aca p=%d nsyn %f %f average %f\n",p,print_parameter/Options.runs,nsyn_persite[p]/Options.runs,average);*/


       aa=*(P+p);
       while(aa==*(P+p)) {
              *(P+p)=MakeMut(pmatrix+((*(P+p))*AA));
      }
      mut++;
  }

  if (aa!=*(P+p)) *nsyn=1;
  position=p;

#if defined(DEBUG)
     printf("pos=%d antes=%d despu=%d eijmat= %f ",p+1,aa+1,*(P+p)+1,eijmat[aa][*(P+p)]);
     position=p;
     printf("pos=%d  ",position);
#endif

  return(mut);
}

int MakeMut (double *pmatrix)
{
  int i;
  double r;

  r=rndu();
  for (i=0; r>(*pmatrix) && i< AA; i++) pmatrix++;
  return (i);
}

double Calculate_score( int *seqnumT, int *seqnumA)
{
  int j,cont,i;
  double score;
  double energyT, energyA;
  double energyTh, energyAh;
  double energyTsum, energyAsum;
  double energyTsumh, energyAsumh;
  int chainnumber;

  score=energyT=energyA=energyAsum=energyTsum=0.;
  energyTh=energyAh=energyAsumh=energyTsumh=0.;

/*Actualize numprot_total for the calculation of contacts */

      *(numprot_totalA+initial-1+position)= *(seqnumA+position);
      *(numprot_totalT+initial-1+position)= *(seqnumT+position);

  cont=0;
  while(cmatvector[position][cont]) {
      *(numprot_totalA+cmatvector[position][cont])= *(seqnumA+cmatvector[position][cont]-initial+1);
      *(numprot_totalT+cmatvector[position][cont])= *(seqnumT+cmatvector[position][cont]-initial+1);
      cont++;
  }

  if(Options.homo){
       chainnumber=total_length/protsize;
       for(i=0;i<chainnumber;i++){
                   cont=0;
                   while(cmatQvector[position][cont]) {
                        *(numprot_totalA+(i*protsize)+cmatQvector[position][cont])= *(seqnumA+cmatQvector[position][cont]-initial+1);
                        *(numprot_totalT+(i*protsize)+cmatQvector[position][cont])= *(seqnumT+cmatQvector[position][cont]-initial+1);
                        cont++;
                   }
       }
       for(i=0;i<chainnumber;i++){
                   cont=0;
                   while(cmatvector[position][cont]) {
                        *(numprot_totalA+(i*protsize)+cmatvector[position][cont])= *(seqnumA+cmatvector[position][cont]-initial+1);
                        *(numprot_totalT+(i*protsize)+cmatvector[position][cont])= *(seqnumT+cmatvector[position][cont]-initial+1);
                        cont++;
                   }
       }
   }

/*        printf("position %d\n ", position);
    for(j=0;j<protsize;j++) 
        printf("%d-%d ", j,*(seqnumA+j));
        printf("\n");
    for(j=0;j<protsize;j++) 
        printf("%d-%d ", j,*(seqnumT+j));
        printf("\n");
    for(j=0;j<total_length;j++) 
        printf("%d-%d ", j,*(numprot_totalA+j));
        printf("\n");
    for(j=0;j<total_length;j++) 
        printf("%d-%d ", j,*(numprot_totalT+j));
        printf("\n");
    for(j=0;j<MAXCONTACTS;j++) printf("%d ", cmatvector[position][j]);
        printf("\n");*/


  if(Options.multimer){
    cont=0;
    while(cmatQvector[position][cont]) {
      energyAh =  eijmat[*(seqnumA+position)][*(numprot_totalA+cmatQvector[position][cont])];
      energyTh =  eijmat[*(seqnumT+position)][*(numprot_totalT+cmatQvector[position][cont])];
 
      energyAsumh+= energyAh;
      energyTsumh+= energyTh;
      cont++;
    }
    cont=0;
    while(cmatvector[position][cont]) {
      energyA =  eijmat[*(seqnumA+position)][*(numprot_totalA+cmatvector[position][cont])];
      energyT =  eijmat[*(seqnumT+position)][*(numprot_totalT+cmatvector[position][cont])];

      energyAsum+= energyA;
      energyTsum+= energyT;

      score += ( energyT-energyA)*( energyT-energyA);
      cont++;
     }

    score += ((energyTsumh+energyTsum)-(energyAsumh+energyAsum))*((energyTsumh+energyTsum)-(energyAsumh+energyAsum));

   }
    else { 
            cont=0;
            while(cmatvector[position][cont]) {
                  energyA =  eijmat[*(seqnumA+position)][*(numprot_totalA+cmatvector[position][cont])];
                  energyT =  eijmat[*(seqnumT+position)][*(numprot_totalT+cmatvector[position][cont])];

                 energyAsum+= energyA;
                 energyTsum+= energyT;

                 score += ( energyT-energyA)*( energyT-energyA);
                 cont++;
           }
           score += (energyTsum-energyAsum)*(energyTsum-energyAsum);
  }

   score=sqrt(score);

   return(score);

}


int Accept_trial(double score,float beta)
{
  double prob,r;

  if(score==0) return(1);


  if(Options.acceptmodel) {       /*hard evaluation*/

#if defined(DEBUG)
     printf("score=%f beta=%f ",score, beta);
#endif
                        if(score<beta) return(1);
                        else return(0);
  }
  else{

     prob=(-beta*(score*score))/(1-(exp(beta*score*score)));    /*soft evaluation*/

     /*Evaluate if trial is accepted or not*/

     r=rndu();

#if defined(DEBUG)
     printf("prob=%f r=%f score=%f beta=%f ",prob, r, score,beta);
#endif
     if(r>prob) return(0);    
     else return(1);
  }
}

char run(int *Aj,int *Tj,int *mut)
{

  int i,j;
  int *A,*T,*P;

  A=Aj;
  T=Tj;

  *mut=0;
  memcpy(T,A,protsize*(sizeof(int)));         /*Copy accepted to trial and mutate it*/
  *mut=MutateSeq(T,pmatrix,&nsyn);

#if defined(DEBUG_SEQ)
  for(i=0;i<protsize;i++) printf("%d %d.%d ",i+1,*(A+i),*(T+i));
  printf("\n");
#endif

  if(*mut==0) {
    Accumul_sust(A,A,Nsus_pos);
    return('0');
  }
  else {
    

    score=Calculate_score(T, A);

    if(score==0) {
                   Accumul_sust(A,T,Nsus_pos);   
                   memcpy(A,T,protsize*(sizeof(int)));    /*if trial is Accepted copy T -> A*/
                   return('A');  /* in screening mode all substitution are nonsynonimous, score could be 0 for example due to contact number is 0.THis site compositon probably varies randomly*/
    }
    else{     /* positive and negative scores*/

      if(Accept_trial(score,beta)){
        Accumul_sust(A,T,Nsus_pos);   
        memcpy(A,T,protsize*(sizeof(int)));    /*if trial is Accepted copy T -> A*/
        return('A');
      }	
      else{
           Accumul_sust(A,A,Nsus_pos);  /* trial Not accepted */
           return('N');
      }
    }
  }
}

void Memory(void)
{

 if((qmatrix=(double *)malloc((AA*AA)*sizeof(double)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((pmatrix=(double *)malloc((AA*AA)*sizeof(double)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((numprot=(int *)malloc(protsize*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((numprot_total=(int *)malloc(total_length*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((numprot_totalA=(int *)malloc(total_length*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((numprot_totalT=(int *)malloc(total_length*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((A=(int *)malloc(Options.runs*protsize*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((T=(int *)malloc(Options.runs*protsize*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((freq=(double *)malloc(AA*sizeof(double)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
  if((Nsus_pos=(int *)malloc(protsize*AA*AA*sizeof(int)))==NULL){
    fprintf(stderr,"Memory allocation failed!!\n");
    exit(1);
  }
                  

}

int Init_vector(int *vector,int vecsize)

{
  int i;

  for(i=0;i<vecsize;i++) *(vector+i)=0;
  return(0);
}

int Accumul_sust(int *seqA,int *seqT,int *Nsus_pos)
{

  int i;

    for(i=0;i<protsize;i++)
      (*(Nsus_pos+(i*AA*AA)+ (*(seqA+i))*AA + *(seqT+i)))++;

  return(0);
}


int Count(char ch,double *nonaccepted, double *accepted, double *nothing,int *mut,int *nsyn,double *N_nsyn,double * N_syn)
{


  switch(ch){
  case 'A': *accepted= *accepted + *mut; break;
  case 'N': *nonaccepted=*nonaccepted + *mut; break;
  case '0': *nothing++; break;
  }

  if(ch=='A'){
                  if(*nsyn!=1) {
                                 (*N_syn)++;
                                  syn_persite[position]++;
                  }
                       else    {
                              (*N_nsyn)++;
                              nsyn_persite[position]++;
                    }
        }
  if(ch=='N')     {
                (*N_syn)++;
                syn_persite[position]++;
  }

  return(0);
}

