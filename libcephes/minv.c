/*							minv.c
 *
 *	Matrix inversion
 *
 *
 *
 * SYNOPSIS:
 *
 * int n, errcod;
 * double A[n*n], X[n*n];
 * double B[n];
 * int IPS[n];
 * int cfs_minv();
 *
 * errcod = cfs_minv( A, X, n, B, IPS );
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the inverse of the n by n matrix A.  The result goes
 * to X.   B and IPS are scratch pad arrays of length n.
 * The contents of matrix A are destroyed.
 *
 * The routine returns nonzero on error; error messages are printed
 * by subroutine simq().
 *
 */

int cfs_minv( A, X, n, B, IPS )
double A[], X[];
int n;
double B[];
int IPS[];
{
double *pX;
int i, k;
extern int cfs_simq(double *, double *, double *, int, int, int *);
extern void cfs_mtransp(int, double *, double *);

for( i=1; i<n; i++ )
	B[i] = 0.0;
B[0] = 1.0;
/* Reduce the matrix and solve for first right hand side vector */
pX = X;
k = cfs_simq( A, B, pX, n, 1, IPS );
if( k )
	return(-1);
/* Solve for the remaining right hand side vectors */
for( i=1; i<n; i++ )
	{
	B[i-1] = 0.0;
	B[i] = 1.0;
	pX += n;
	k = cfs_simq( A, B, pX, n, -1, IPS );
	if( k )
		return(-1);
	}
/* Transpose the array of solution vectors */
cfs_mtransp( n, X, X );
return(0);
}

