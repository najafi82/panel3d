#define MAX(a,b) ((a > b) ? (a) : (b))
#define MIN(a,b) ((a < b) ? (a) : (b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) ((a == 0.0) ? 0.0 : a*a)

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)calloc(nh-nl+2,sizeof(double));
	if (!v)
    {
        puts("allocation failure in dvector()");
        exit(-1);
	}
	return v-nl+1;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) calloc(nrow+1,sizeof(double*));
	if (!m)
    {
        puts("allocation failure 1 in matrix()");
        exit(-1);
	}
	m += 1;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) calloc(nrow*ncol+1,sizeof(double));
	if (!m[nrl])
	{
        puts("allocation failure 2 in matrix()");
        exit(-1);
	}
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++)
        m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free(v+nl-1);
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
	free (m[nrl]+ncl-1);
	free (m+nrl-1);
}

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmp(double **a, int m, int n, double w[], double **v)
{
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=dvector(1,n);
	g=0.0;
	scale=0.0;
	anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=0.0;
		s=0.0;
		scale=0.0;
		if (i <= m)
		{
			for (k=i;k<=m;k++)
                scale += fabs(a[k][i]);
			if (scale)
			{
				for (k=i;k<=m;k++)
				{
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++)
				{
				    s=0.0;
					for (k=i;k<=m;k++)
                        s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++)
                        a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++)
                    a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=0.0;
		s=0.0;
		scale=0.0;
		if (i <= m && i != n)
        {
			for (k=l;k<=n;k++)
                scale += fabs(a[i][k]);
			if (scale)
			{
				for (k=l;k<=n;k++)
				{
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++)
                    rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++)
				{
					for (s=0.0,k=l;k<=n;k++)
                        s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++)
                        a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++)
                    a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--)
    {
		if (i < n)
		{
			if (g)
			{
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++)
				{
					for (s=0.0,k=l;k<=n;k++)
                        s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++)
                        v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++)
                v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n);i>=1;i--)
	{
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++)
            a[i][j]=0.0;
		if (g)
		{
			g=1.0/g;
			for (j=l;j<=n;j++)
			{
				for (s=0.0,k=l;k<=m;k++)
                    s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++)
                    a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++)
                a[j][i] *= g;
		}
        else
            for (j=i;j<=m;j++)
                a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--)
    {
		for (its=1;its<=300;its++)
		{
			flag=1;
			for (l=k;l>=1;l--)
			{
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm)
				{
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm)
                    break;
			}
			if (flag)
			{
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++)
				{
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm)
                        break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++)
					{
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 300)
			{
                puts("no convergence in 300 svdcmp iterations");
                exit(-1);
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=1.0;
            s=1.0;
			for (j=l;j<=nm;j++)
            {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++)
				{
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z)
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++)
				{
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s,*tmp;

	tmp=dvector(1,n);
	for (j=1;j<=n;j++)
    {
		s=0.0;
		if (w[j])
		{
			for (i=1;i<=m;i++)
                s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++)
	{
		s=0.0;
		for (jj=1;jj<=n;jj++)
            s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_dvector(tmp,1,n);
}

int SOR(struct MatrixCoefficient *Matrix,double **Rhs,double **X,int NoRow,double Eps,int MaxIter,double Omega,double *Error)
{
    int i,k=0,j;
    double S,error=0.0,error1=0.0,zigma;
    do
    {
        k++;
        error=0.0;
        for(i=0;i<NoRow;i++)
        {
            S=0.0;
            for (j=0;j<NoRow;j++)
                if(i!=j)
                    S+=Matrix->Elem[i][j]* (*X)[j];

            zigma=((*Rhs)[i]-S)/Matrix->Elem[i][i];
            error1=zigma-(*X)[i];
            if (fabs(error1)>error)
                error=fabs(error1);
            (*X)[i]=(1.0-Omega)* (*X)[i] + Omega*zigma;
        }
    }while((error>Eps)&&(k<MaxIter));
    *Error=error;
    return (k);
}

