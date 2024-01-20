

#include "mex.h"
#include "math.h"


void inerciatoep_c(double *A, double *G,  double *H, double *q, mwSize n, mwSize alfa)
{ //A is n*n, G is n * alfa
    mwSize i,j,m,k;
    double *f, *fn, *y,aux, *vaux1, *vaux2, *w, *paux;
    //
    f=(double*)malloc(n*alfa*sizeof(double));
    fn=(double*)malloc(n*alfa*sizeof(double));
    y=(double*)malloc(n*sizeof(double));
    vaux1=(double*)malloc(alfa*sizeof(double));
    vaux2=(double*)malloc(n*sizeof(double));
    w=(double*)malloc(n*n*sizeof(double));
    q[0]=A[0];
    w[0]=A[n]/q[0]; //A[n]=A(1,2)
   
    
    for (i=0;i<alfa;i++)
     {
        f[0+i*n]=G[0+i*n]/q[0];
      }
    for (m=1; m<n-1 ;m++)
    {
        //memset(fn, 0, 0*sizeof(double)); //fn=zeros(m,alfa);
        q[m]=A[m+n*m];
        for (j=0;j<m;j++)
            q[m]-=A[j+n*m]*w[j+n*(m-1)];  
        for (j=0;j<m;j++)
            y[j]=w[j+n*(m-1)];
        y[m]=-1;
        for (i=0;i<alfa;i++)
           {aux=0;
            for (k=0;k<m;k++)
               aux+=A[k+n*m]*f[k+n*i];
            // aux = vmm1'*f(1:m-1,i)
              
             for(j=0;j<m;j++)
             {             
               fn[j+i*n]= f[j+i*n]-(G[m+n*i]-aux)*y[j]/q[m];
             }
            fn[m+i*n]= 0-(G[m+n*i]-aux)*y[m]/q[m];
           }
        for (i=0;i<alfa;i++)
           {vaux1[i]=0;
           for (j=0;j<=m;j++)
               vaux1[i]+=H[j+i*n]*y[j];
           }
        for (i=0;i<=m;i++)
           {vaux2[i]=0;
            for (j=0;j<alfa;j++)
                vaux2[i]+=fn[i+j*n]*vaux1[j];
            }
        w[0+n*m]=0-vaux2[0];
        for(j=1;j<=m;j++)
            w[j+n*m]=w[j-1+n*(m-1)]-vaux2[j];
        paux=f;
        f=fn;
        fn=paux;
        
    }
     q[n-1]=A[n-1+n*(n-1)];
     for (j=0;j<n-1;j++)
            q[m]-=A[j+n*(n-1)]*w[j+n*(n-2)];  
    free(f);
    free(fn);
    free(y);
    free(vaux1);
    free(vaux2);
    free(w);
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *A, *G, *H, *q, *w;              
    double *D, *P, *F;
    size_t nA,n_alfa;            

    int i;
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","3 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
    }
   
    
   A = mxGetPr(prhs[0]);
   G = mxGetPr(prhs[1]);
   H = mxGetPr(prhs[2]);
   
   nA = mxGetM(prhs[0]);
   n_alfa=mxGetN(prhs[1]);
        
  
    plhs[0] = mxCreateDoubleMatrix((mwSize)nA,(mwSize)1,mxREAL); //D
    q = mxGetPr(plhs[0]);
     
   inerciatoep_c(A, G,H, q, nA, n_alfa) ;   
  
}


