#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
void endrun();
struct io_header
{
  int npart[6];                        /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
                                            stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                        /*!< flags whether the simulation was including star formation */
  int flag_feedback;                   /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];          /*!< total number of particles of each type in this snapshot. This can be
                                            different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;                    /*!< flags whether cooling was included  */
  int num_files;                       /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;                 /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                     /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
  char fill[60];                       /*!< fills to 256 Bytes */
}
 header;

FILE *snapfile, *fp;
char snapfname[20],SN[5];
int check2,l,len,result=0;
float check;

int main(int argc, char **argv){
  
  if(argc !=3){
    printf("./a.out number_of_snapshots snapshot_\n");
    return(-1);
  }//if

// gets(argv[1]);

len = strlen(argv[1]);
float time=0, Amax=0;
for(l=0; l<len; l++)
  {
 result = result * 10 + ( argv[1][l] - '0' );
  }
//printf("%d \n \n",result);


for(l=0;l<result;l++)

{ 
   sprintf(SN,"%03d",l);
     // printf("%s",SN);
   sprintf(snapfname,"%s",argv[2]);
    
   strcat(snapfname,SN);
     
 //printf("%s\n\n",snapfname);

 // printf("Opening %s \n",snapfname);   
  snapfile = fopen(snapfname,"r");

  int blockstart,blockend;
  fread(&blockstart,sizeof(int),1,snapfile);
  fread(&header,sizeof(struct io_header),1,snapfile); 
  fread(&blockend,sizeof(int),1,snapfile);
  if(blockstart!=blockend)
    endrun();

  int i,j;
  int ndim =3 ;
  int nptype =6;
//  for(i=0;i<nptype;i++)
 //  printf("npart[%d]=%d\n",i,header.npart[i]);
 // printf("BoxSize=%g Om0=%g Oml=%g h=%g \n",header.BoxSize,header.Omega0,header.OmegaLambda,header.HubbleParam);
   
 // printf("blockstart=%d blockend=%d\n",blockstart,blockend);


  float *pos;
  int np=0;  
  for(i=0;i<nptype;i++)
  np+=header.npart[i];
  pos = (float *)malloc(sizeof(float)*np*ndim);
  fread(&blockstart,sizeof(int),1,snapfile);
  int blockcheck=(int)sizeof(float)*np*ndim;
  fread(pos,sizeof(float),ndim*np,snapfile);
  fread(&blockend,sizeof(int),1,snapfile);

 // printf("blockstart=%d blockcheck= %d blockend= %d \n\n",blockstart,blockcheck,blockend);
  //return(-1); 

  float *vel;
  vel= (float *)malloc(sizeof(float)*np*ndim);
  fread(&blockstart,sizeof(int),1,snapfile);
  fread(vel,sizeof(float),ndim*np,snapfile);
  fread(&blockend,sizeof(int),1,snapfile);
 // printf("blockstart=%d blockend= %d \n\n",blockstart,blockend);
  
// for(i=0;i<np;i++) 
//  printf("%f %f  %f \n ",vel[0+3*i],vel[1+3*i],vel[2+3*i]);
//  printf("blockstart=%d blockcheck= %d blockend= %d \n",blockstart,blockcheck,blockend);
  

 unsigned int *ID;
 ID=(unsigned int *)malloc(sizeof(unsigned int)*np);
 fread(&blockstart,sizeof(int),1,snapfile);
 fread(ID,sizeof(unsigned int),np,snapfile);
 fread(&blockend,sizeof(int),1,snapfile);

//caculating com position and velocity of particles
float Vcm[3]={0},Pcm[3]={0},Vcmh[3]={0},Pcmh[3]={0},Vcmd[3]={0},Pcmd[3]={0};
for(i=0;i<2000000;i++)
{
 if(i<1000000)
{  
   Vcmh[0]+=vel[0+i*3];
   Vcmh[1]+=vel[1+i*3];
   Vcmh[2]+=vel[2+i*3];
   Pcmh[0]+=pos[0+i*3];
   Pcmh[1]+=pos[1+i*3];
   Pcmh[2]+=pos[2+i*3];
}
else
  {
   Vcmd[0]+=vel[0+i*3];
   Vcmd[1]+=vel[1+i*3];
   Vcmd[2]+=vel[2+i*3];
   Pcmd[0]+=pos[0+i*3];
   Pcmd[1]+=pos[1+i*3];
   Pcmd[2]+=pos[2+i*3];
  }
}

Vcmd[0]=Vcmd[0]/1000000;            // Center of mass and it's velocity for Disk particles
Vcmd[1]=Vcmd[1]/1000000;
Vcmd[2]=Vcmd[2]/1000000;
Pcmd[0]=Pcmd[0]/1000000;
Pcmd[1]=Pcmd[1]/1000000;
Pcmd[2]=Pcmd[2]/1000000;

Vcmh[0]=Vcmh[0]/1000000;             // Center of mass and it's velocity for halo particles       
Vcmh[1]=Vcmh[1]/1000000;
Vcmh[2]=Vcmh[2]/1000000;
Pcmh[0]=Pcmh[0]/1000000;

//printf("%f %f %f \n" , Pcmd[0],Pcmd[1],Pcmd[2]); 
//Disk particles in frame of cm_disk reference

for(i=0;i<2000000;i++)
{
  if(i>999999)
  {  
   pos[0+i*3]=pos[0+i*3]-Pcmd[0];
   pos[1+i*3]=pos[1+i*3]-Pcmd[1];
   pos[2+i*3]=pos[2+i*3]-Pcmd[2];

   vel[0+i*3]=vel[0+i*3]-Vcmd[0];
   vel[1+i*3]=vel[1+i*3]-Vcmd[1];
   vel[2+i*3]=vel[2+i*3]-Vcmd[2];
  }
}



// Calculating bar strength

// 1. Radius of particles
float R[2000000];

for(i=0;i<2000000;i++)
{
   R[i]=sqrt(pos[0+i*3]*pos[0+i*3] +pos[1+i*3]*pos[1+i*3]);
// if (i>999999 & i<1000300)
 //  printf("%f \n", R[i]);
}

int nstar=header.npart[3],m=0;
float a[3]={0.0,0.0,0.0},b[3]={0.0,0.0,0.0},theta=0,A2=0,A2max=0,R_diff,phi_1,phi_2,phi_diff;

if(l==0){
  phi_2=0;
}
for(R_diff=0;R_diff<35;R_diff++){
  for(i=1000000;i<2000000;i++){
 //      if(R[i]>R_diff && R[i]<R_diff+1 && pos[3+3*i]<1.5 && pos[3+3*i]>-1.5)
    if(R[i]>R_diff && R[i]<R_diff+1){
         theta=atan(pos[1+3*i]/pos[0+3*i]);
         a[0]+=1;
         //   a[1]+=cos(theta);
         a[2]+=cos(2*theta);  //Real Part of fourier mode
         //   b[1]+=sin(theta);
         b[2]+=sin(2*theta);  // imaginary part of fourier mode
    }
  }

  if (R_diff==4){
    phi_1=0.5*atan(b[2]/a[2]);
  }
  //A2=(sqrt(a[2]*a[2]+b[2]*b[2])/a[0]);
  //printf("A2=%f \n",A2); 
  if(A2max<A2){ 
     A2max=A2;
     check=R_diff;
  }
}
//printf("check \n");

if(phi_1 >0 && phi_2>0){
  phi_diff=phi_1 - phi_2;
}
else if (phi_1>0 && phi_2<0){
  phi_diff=phi_1-phi_2;
}
else if (phi_1<0 && phi_2 >0){
 phi_diff=phi_1-phi_2+1.57;
}
else if (phi_1 <0 && phi_2<0){ 
 phi_diff=phi_2 - phi_1 ;
}

phi_2=phi_1;

if(phi_diff<0)
 phi_diff=-phi_diff;

printf("\n %d  %f %f \n",l,(float)l*0.0196,phi_diff/0.0196); // pattern speed in unit of km/s /kpc

  free(pos);
  free(vel);
  free(ID);
  fclose(snapfile);
  
}//for_after_main 

//printf ("snapshot=%d Rmax=%f Amax=%f \n", check2,check,Amax);  

return(-1);

}//main() 


void endrun(){
   exit(-1);
}
