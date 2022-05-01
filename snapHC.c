//input: lattice size (in unit cells. Total number of sites is 2*L^2) and file with data
//output: data converted to the hexagonal/honeycomb lattice


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>

int createLattice(int);
void plotEps(int, char*);

int main(int argc, char *argv[])
{
  FILE *fp, *out;
  int x,y,est,lix, LSIZE;
  double c30,s30,xf,yf,lxf ;

  //pass size L and file as command line arguments  
  if(argc!=3)
    {printf("./a.out L file\n"); return 0;}
  
  LSIZE = atoi(argv[1]);
  c30 = cos(M_PI/6.0);
  s30 = sin(M_PI/6.0);  
  fp=fopen(argv[2],"r");
  
  out=fopen("convertedData.cfg", "w");

  createLattice(LSIZE);
  
  int est2;
  //the input file has 4 columns. If the unit cell is i=0...L^2 and the basis atom in the unit cell is j=0 or 1, then
  // y = i/L
  // x = 2*(i%L) + j - i/L
  // est = 0 if empy, 1 if occupied by particle

  //this loop will convert the data to the proper coordinates
  while((fscanf(fp,"%d %d %d %d",&lix,&y,&x,&est)==4))
    {
      x++;
      xf=x*c30;
      if(y%2)
	{
	  if(x%2)	    
	    yf = y*(1+s30);	    
	  else
	    yf = (y*3.+1.)/2.;	    	  
	}     
      else
	{	  
	  if(x%2)
	    yf = 1.5*y+0.5;
	    
	  else	    
	    yf = 1.5*y;
	}    
      fprintf(out,"%f %f %d\n",xf,yf,est);
    }
  fclose(fp);
  fclose(out);
  
  //plotEps(LSIZE, argv[2]);
  return 0;
}
  

//this function will create a file with the points to draw the lattice lines.
//The points are constructed in a way suitable to gnuplot
#define Z 1.0
int createLattice(int LSIZE)
{
  int i,j;
  float x0,y0,x1,y1;
  float PI6,PI3;
  FILE *arq;
  char nFile[20];
  sprintf(nFile, "./%d.hc", LSIZE);
  arq=fopen(nFile, "w");

  int LL=3*LSIZE;
  
  PI6 = M_PI/6.0;
  PI3 = M_PI/3.0;
  
  x0 = 0.0;
  y0 = 0.0;

  for(j=0;j<(LL-1);j=j+2) {
    for(i=0;i<LL;i++) {
      x1 = x0 + i*Z*cos(PI6);
      y1 = y0 + (i%2)*Z*sin(PI6);
      fprintf(arq, "%f %f\n",x1,y1);
    }
    fprintf(arq, "\n\n");
    for(i=1;i<LL;i=i+2) {
      x1 = x0 + i*Z*cos(PI6);
      y1 = y0 + (i%2)*Z*sin(PI6);
      fprintf(arq, "%f %f\n",x1,y1);
      y1 = y0 + (i%2)*Z*sin(PI6) + Z;
      fprintf(arq, "%f %f\n",x1,y1);
      fprintf(arq, "\n\n");
    }
    fprintf(arq, "\n\n");
    
    y0 = y0 + (2*sin(PI6)+1)*Z;
    for(i=0;i<LL;i++) {
      x1 = x0 + i*Z*cos(PI6);
      y1 = y0 - (i%2)*Z*sin(PI6);
      fprintf(arq, "%f %f\n",x1,y1);
    }
    fprintf(arq, "\n\n");
    for(i=0;i<LL;i=i+2) {
      x1 = x0 + i*Z*cos(PI6);
      y1 = y0 + (i%2)*Z*sin(PI6);
      fprintf(arq, "%f %f\n",x1,y1);
      y1 = y0 + (i%2)*Z*sin(PI6) + Z;
      fprintf(arq, "%f %f\n",x1,y1);
      fprintf(arq, "\n\n");
    }
    y0 = y0 + Z;    
  }
  fclose(arq);
  return 0;  
}

//calls the gnuplot script (can also be done manually)
void plotEps(int LSIZE, char* arq)
{
  char com[50];
  sprintf(com,"gnuplot -e \"L=%d; file='%s'\" plotSnap.gp", LSIZE, arq);
  system(com); 
  system("rm convertedData.cfg && rm *.hc");
}
