#include "stdlib.h"

/****************************************/
  void OrVector(int *n, int *x, int *y)
    {
      int i;
      for(i=0; i<(*n); i++)
    {
          if(*(x+i) ==1)
        *(y+i)=1;
          }
    }

/*****************************************/
void subtree(int *n, int *merge, int *node, int *pknode, int *pktree)
{
  int i,j,k,temp,length,temp1,temp2, nn;
  int *letree, *ritree;
  int *lenode, *rinode;

  temp=*node;
  length=*n;  
  
  letree=malloc((length+1)*sizeof(int));
  ritree=malloc((length+1)*sizeof(int));
  lenode=malloc(length*sizeof(int));
  rinode=malloc(length*sizeof(int));
  
  nn=length+1;
  i=temp-1;
  temp1=*(merge+i);
  j=length+i;
  temp2=*(merge+j);

  for(i=0;i<nn;i++)
    {
      letree[i]=0;
      ritree[i]=0;
    }
  for(j=0;j<length;j++)
    {
      lenode[j]=0;
      rinode[j]=0;
    }
  
  if(temp1<0) 
    pktree[-temp1-1]=1;
  else
    {
    subtree(n, merge, &temp1, lenode, letree);
    OrVector(n,lenode, pknode);
    OrVector(&nn, letree, pktree);
    }

  if(temp2<0) 
    pktree[-temp2-1]=1;
  else
    {
    subtree(n, merge, &temp2, rinode, ritree);
    OrVector(n,rinode, pknode);
    OrVector(&nn, ritree, pktree);
    }

  *(pknode+temp-1)=1;
  
  free(letree);
  free(ritree);
  free(lenode);
  free(rinode);
}
     
/*************************************************/
void setzero(int *n, int *x)
{
  int i;
  for(i=0; i<(*n); i++)
     {
       *(x+i)=0;
     }
}

void maxnonzero(int *n, int *x, int *place)
{
  int i;
  *place=0;
  for(i=(*n)-1; i>=0; i--)
    {
      if(*(x+i)>0)
    {
      *place=i+1;
          break;
    }
    }
}

/**********************************************************
For a tree with n+1 leaves (n nodes), find the corresponding leaves to 
a certain subset of selected nodes.
senode:    vector of length n;   1--selected, 0-notselected
seleaf:    vector of length n+1;    index of clusters 
***********************************************************/

void selectclust(int *n, int *merge, int *senode, int *seleaf)
{
  int i,j,k,temp,nn;
  int *sn, *tleaf, *tnode;
  
  nn=(*n)+1;
  sn=malloc((nn-1)*sizeof(int));
  tleaf=malloc(nn*sizeof(int));
  tnode=malloc(nn*sizeof(int));
  
  setzero(&nn, tleaf);
  setzero(n, tnode);

  for(i=0;i<(*n);i++)
    {
      *(sn+i)=*(senode+i);
    }

  maxnonzero(n, sn, &temp);
  k=1;
  while(temp>0)
    {
      subtree(n, merge, &temp, tnode, tleaf);
      for(i=0;i<((*n)+1);i++)
    {
      if(*(tleaf+i)==1)
        *(seleaf+i)=k;
    }
      for(i=0;i<(*n);i++)
        {
          if(*(tnode+i)==1)
        *(sn+i)=0;
    }
      k++;
      setzero(&nn, tleaf);
      setzero(n, tnode);
      maxnonzero(n, sn, &temp);
    }
   free(sn);
   free(tleaf);
   free(tnode); 
}
  


/********************************************************
clactree:   big vector M by 7
            1: size
            2: height
            3: mean
            4,5: merge
            6: chrom
            7: arm
ch, nuc:    vector of length N
chset:      vector of length k
senode:     vector of length M
seleaf:     vector of length N
*********************************************************/
void selectleaf(int *M, int *clactree, int *N, int *ch, int *k, int *chset, int *nuc, int *centro, int *senode, int *seleaf, int *cc)
{
  int i,j,nn,kk,mm;
  int *mergetemp;
  int *senodetemp, *seleaftemp;
  int curch, curcen, arm;
  int curpoint, count, a,b,c,d, ntemp;
  int begin, beginlabel, end;
  

  kk=*k;
  mm=*M;
  nn=*N;
  
  senodetemp=malloc(mm*sizeof(int));
  seleaftemp=malloc(mm*sizeof(int));
  mergetemp=malloc(mm*2*sizeof(int));
  
  curpoint=0;
  for(i=0; i<kk; i++)
      {
        curch=*(chset+i);
    curcen=*(centro+curch-1);
        /* arm=1, left; arm=2, right */
        for(arm=1;arm<3;arm++)
      {
      begin=0;
          beginlabel=0;
          end=0;
          for(j=0; j<mm; j++)
        {
              if(*(clactree+mm*2+j) ==curch & *(clactree+mm*3+j)==arm)
                {
                  if(beginlabel==0)
            {
                    begin=j;
                    beginlabel=1;
                    }
          end=j;
        }
        }
      /******** take out the tree from "begin" to "end" ****/
          if(beginlabel==1)
          {  
             ntemp=end-begin+1;
         for(a=0;a<ntemp;a++)
              {
                  mergetemp[a]=*(clactree+a+begin);
                  mergetemp[ntemp+a]=*(clactree+mm+a+begin);
                  senodetemp[a]=*(senode+a+begin);
                  seleaftemp[a]=0;
              }
             seleaftemp[ntemp]=0;
             selectclust(&ntemp, mergetemp, senodetemp, seleaftemp); 
         for(a=0;a<=ntemp;a++)
                   {
             *(seleaf+curpoint+a)=*(seleaftemp+a);
           }
             curpoint=curpoint+ntemp+1;
           }
      else
        {
          count=0;
             for(b=0; b<nn; b++)
           {
                if(*(ch+b) ==curch & ((*(nuc+b)<curcen & arm==1) | (*(nuc+b)>curcen & arm==2)))
                  count=count+1;
        }
             if(count>0)
               {
                 for(b=0; b<count;b++)
                   {
             *(seleaf+curpoint+b)=0;
                   }
                 curpoint=curpoint+count;
               }
        }
      }
      }
  *cc=curpoint;
  
  free(senodetemp);
  free(seleaftemp);
  free(mergetemp);
}              
                  
/***********************************************************************
 ***********************************************************************
calculate region mean value
sample:     vector of length n
result:     vector of length n

clustindex: vector of length n

ch:         vector of length n
chset:      vector of length k
nucpos:     vector of length n
************************************************************************/
void makemean(int *n, double *sample, int *clustindex, int *ch, int *chset, int *k, int *nucpos, int *centro, double *result)
{
  int i,j,l,m, nn, kk;
  int curch, curcen;  
  int begin, beginlabel, end;
  int totalClust, arm;
  int count;
  double sumtemp;

  nn=*n;
  kk=*k;
  for(i=0; i<kk; i++)
    {
      curch=*(chset+i);
      curcen=*(centro+curch-1);
     for(arm=1;arm<3;arm++)
       { 
      begin=0;
      beginlabel=0;
      end=0;
      totalClust=0;
      for(j=0;j<nn;j++)
    {
          
          if(*(ch+j)==curch & ( (*(nucpos+j)<curcen & arm==1) | (*(nucpos+j)>curcen & arm==2)))
            {
              if(beginlabel==0)
        {
          begin=j;
                  beginlabel=1;
                }
          end=j;
              if(*(clustindex+j)>totalClust)
        {
                  totalClust=*(clustindex+j);
            }
        }
    }
      /********** calculate mean for each cluster ***/
      if(totalClust>0)
    {
          for(l=1; l<=totalClust; l++)
            {
              sumtemp=0;
              count=0;
              for(m=0;m<(end-begin+1);m++)
                {
                  if(*(clustindex+begin+m)==l)
                    {
                      sumtemp=sumtemp+ *(sample+begin+m);
                      count=count+1;
                    }
        }
              if(count>0)
        {
                  sumtemp=sumtemp/count;
                }
              for(m=0;m<(end-begin+1);m++)
                {
                  if(*(clustindex+begin+m)==l)
                    {
                      *(result+begin+m)=sumtemp;
                    }
        }
        }
    }  
       }
    }
}
