#include "stdio.h"
#include "stdlib.h"

double myabs(double x)
{
  if(x<0)
    {
    return(-x);
    }
  else
    {
    return(x);
    }
}

void testabs(double *x, double *y)
{
  double temp;
  temp=*x;
  *y=myabs(temp);
}

/***************************/
double relativediff(double x, double y)
{  
  double temp, z; 
  temp=myabs(x)+myabs(y)+myabs(x+y)+0.001;  
  z=myabs(x-y)/temp;
  return(z);
}

void testrd(double *x, double *y, double *z)
{
  double temp1, temp2;
  temp1=*x;
  temp2=*y;
  *z=relativediff(temp1, temp2);
  printf("test");
}

/**************************/
int whichismin(int n, double *x)
{
  int i;
  int result;
  double cur;

  result=0;
  for(i=1; i<n; i++)
    {
      cur=*(x+i);
      if(cur < *(x+result)) 
        {
    result=i;
    }
    }
  return (result);
}

void testwim(int *n, double *x, int *p)
{
  int nn;
  nn=*n;
  *p=whichismin(nn, x);
}

/***************************/
void shiftforwardint( int n, int *x, int position)
{
  int i;
  for(i=position; i<(n-1); i++)
    {
    *(x+i)=*(x+i+1);
    }
}

void testsfint(int *n, int *x, int *p)
{
  int nn, pp;
  nn=*n;
  pp=*p;
  shiftforwardint(nn, x, pp);
}  


void shiftforwarddouble( int n, double *x, int position)
{
  int i;
  for(i=position; i<(n-1); i++)
    {
    *(x+i)=*(x+i+1);
    }
}
void testsfdouble(int *n, double *x, int *p)
{
  int nn, pp;
  nn=*n;
  pp=*p;
  shiftforwarddouble(nn, x, pp);
}  


/*****************************************************************
Build one tree for one chrom arm
n:        length of x
x:        array measurement

result:   vector of length total*7
begin:    the beginning row after writing result as a (total by 7) matrix
total:    the total number of rows for result after converting it into a matrix of 7 column

curch:    current chromosome number
arm:      current arm index 1--left, 2--right
******************************************************************/

void hclust_inside(int *n, double *x, double *result, int begin, int total, int curch, int arm)
{
  int i, o, nn;
  int length; 
  int *val1;
  double *val2, leftsum, rightsum;
  double leftsize, rightsize;


  length=*n;
  /* nn=length-1;*/
  nn=total;

  val1=malloc(length*sizeof(int));
  val2=malloc(length*sizeof(double)); 

  /* val1: -1 -- -n =length */
  for(i=0;i<(length-1);i++)
    {
      val1[i]=-i-1;
      val2[i]=relativediff(*(x+i), *(x+i+1)) ;
    }

  val1[length-1]= -length;
  val2[length-1]= 100000;

  for(i=0;i<(length-1);i++)
    {
      /* o <- which.is.min(val[, 2]) */
      o=whichismin(length - i,val2);

      /*  height[ii] <- val[o, 2] */   
      *(result+nn*1+i+begin)=val2[o];

      /* merge[ii,  ] <- c(val[o, 1], val[o + 1, 1]) */
      *(result+nn*3+i+begin)=val1[o];
      *(result+nn*4+i+begin)=val1[o+1];      

      *(result+nn*5+i+begin)=curch;
      *(result+nn*6+i+begin)=arm;

      if(val1[o]<0)
    {
          leftsize=1;
      /*  leftsum<-x[-val[o,1]] */
          leftsum=*(x-val1[o]-1);
    }
      if(val1[o]>0)
        {
          /*leftsize<-size[val[o,1]]*/
          leftsize=*(result+val1[o]-1 + begin);
      /*leftsum<-sumvalue[val[o,1]]*/
          leftsum=*(result+nn*2+val1[o]-1 +begin);
    }
      if(val1[o+1]<0)
        {
          rightsize=1;
          rightsum=*(x-val1[o+1]-1);
    }
      if(val1[o+1]>0)
        {
          /* rightsize<-size[val[o+1,1]] */
          rightsize=*(result+ begin+ val1[o+1]-1);
          /* rightsum<-sumvalue[val[o+1,1]] */
          rightsum=*(result+ nn*2+ begin+ val1[o+1]-1);
        }   
      /* size[ii]<-leftsize+rightsize           */
      *(result+i+begin)=leftsize+rightsize;
      *(result+nn*2+i+begin)=(leftsum+rightsum);
      
      val2[o]=val2[o+1];  
      val1[o]=i+1;
      shiftforwardint(length-i, val1, o+1);
      shiftforwarddouble(length-i, val2, o+1);
    }

  for(i=0; i<(length-1); i++)
    {
    *(result+nn*2+i+begin)=(*(result+nn*2+i+begin))/(*(result+i+begin));
    }
    
    free(val1);
    free(val2);
}


/***********************************
Build the trees for one array
N:       length of the array       
array:   log2 ratio of the array   ;vector of length N

ch:      chromosome number         ;vector of length N
k:       length of chset
chset:   different chrom number    ;vector of length k

nuc:     nucletide position        ;vector of length N
centro:  centromere position       ;vector of length 24/23

result:  returned trees            ;vector of length 7*N
M:       length of effective rows in result
************************************/
void clacarray(int *N, double *array, int *ch, int *k, int *chset, int *nuc, int *centro, double *result, int *M)
{
  int i, j,arm, count;
  int nn, kk, curch, curcen;
  double  *curArm;

  nn=*N;  
  curArm=malloc(nn*sizeof(int));
  
  kk=*k;
  *M=0;
  /* build tree for each chromosome*/
  for(i=0; i<kk; i++)
    {
      curch=*(chset+i);
      curcen=*(centro+curch-1);
      /* arm=1, left; arm=2, right */
      arm=1;
         count=0;
         for(j=0;j<nn; j++)
     {
           if(*(ch+j)==curch & *(nuc+j)<curcen )
             {
           count=count+1;
               curArm[count-1]=*(array+j);
             }
     }
         
     if(count>5)
       {/*build tree for this arm */
             hclust_inside(&count, curArm, result, *M, nn, curch,arm);
           *M=*M+count-1;
           }           
     
      arm=2;
         count=0;
         for(j=0;j<nn; j++)
     {
           if(*(ch+j)==curch & *(nuc+j)> curcen )
             {
           count=count+1;
               curArm[count-1]=*(array+j);
             }
     }
         
     if(count>5)
       {/*build tree for this arm */
             hclust_inside(&count, curArm, result, *M, nn, curch,arm);
           *M=*M+count-1;
           }           
    }
    free(curArm);
}
