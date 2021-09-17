/*
 * configurations
 *--partition=physical
 *--cpus-per-task = 8
 *--mem=32G
 *export OMP_PLACES="cores"
 *export OMP_PROC_BIND="true"
 * */
// CPP program to solve the sequence alignment
// problem. Adapted from https://www.geeksforgeeks.org/sequence-alignment-problem/ 
#include <sys/time.h>
#include <string>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <math.h>
#include <numa.h>
#include "sha512.hh"

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

using namespace std;

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap, int *penalties);
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans);

/*
Examples of sha512 which returns a std::string
sw::sha512::calculate("SHA512 of std::string") // hash of a string, or
sw::sha512::file(path) // hash of a file specified by its path, or
sw::sha512::calculate(&data, sizeof(data)) // hash of any block of data
*/

// Return current wallclock time, for performance measurement
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char **argv){
	int misMatchPenalty;
	int gapPenalty;
	int k;
	//omp_set_nested(1);
	//setenv("OMP_PLACES","threads(1)",1);
	//setenv("OMP_PLACES", "{0:2},{2:2},{4:2},{6:2},{32:2},{34:2},{36:2},{38:2},{8:2},{10:2},{12:2},{14:2},{40:2},{42:2},{46:2},{44:2},{48:2},{50:2},{52:2},{54:2},{16:2},{18:2},{20:2},{22:2},{24:2},{26:2},{30:2},{28:2},{56:2},{58:2},{60:2},{62:2}",1);
	//setenv("OMP_PLACES","{0:8,32:8,8:8,40:8},{48:8,16:8,30:8,56:8}",1);
	//setenv("OMP_PROC_BIND","true",1);
	int nplaces = omp_get_num_places();
	printf("places:%d:pus%d\n",nplaces,omp_get_place_num_procs(1));
	printf("OMP_PLACES:%s\n",getenv("OMP_PLACES"));
	printf("nodes:%d core:%d\n",numa_num_configured_nodes(),numa_num_configured_cpus());
	std::cin >> misMatchPenalty;
	std::cin >> gapPenalty;
	std::cin >> k;	
	std::string genes[k];
	for(int i=0;i<k;i++) std::cin >> genes[i];

	int numPairs= k*(k-1)/2;

	int penalties[numPairs];
		
	uint64_t start = GetTimeStamp ();

	// return all the penalties and the hash of all allignments
	std::string alignmentHash = getMinimumPenalties(genes,
		k,misMatchPenalty, gapPenalty,
		penalties);
		
	// print the time taken to do the computation
	printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
		
	// print the alginment hash
	std::cout<<alignmentHash<<std::endl;

	for(int i=0;i<numPairs;i++){
		std::cout<<penalties[i] << " ";
	}
	std::cout << std::endl;
	return 0;
}

int min3(int a, int b, int c) {
	if (a <= b && a <= c) {
		return a;
	} else if (b <= a && b <= c) {
		return b;
	} else {
		return c;
	}
}

// equivalent of  int *dp[width] = new int[height][width]
// but works for width not known at compile time.
// (Delete structure by  delete[] dp[0]; delete[] dp;)
int **new2d (int width, int height)
{
	int **dp = new int *[width];
	size_t size = width;
	size *= height;
	int *dp0 = new int [size];
	if (!dp || !dp0)
	{
	    std::cerr << "getMinimumPenalty: new failed" << std::endl;
	    exit(1);
	}
	dp[0] = dp0;
	for (int i = 1; i < width; i++)
	    dp[i] = dp[i-1] + height;

	return dp;
}

std::string getMinimumPenalties(std::string *genes, int k, int pxy, int pgap,
	int *penalties)
{
	int probNum=0;
	int x=0;
	std::string alignmentHash="";
	int nmax = k*(k+1)/2;
	std::string allHash[nmax];
	//#pragma omp parallel for schedule(dynamic) num_threads(1) reduction(+:probNum) shared(allHash,nmax) private(x) proc_bind(close)
	for (x=0; x<nmax; x++){
			//printf("outer thread:%d\n",omp_get_thread_num());
			//printf("%d/%d\n",x,nmax);			
			int i = (-1 +sqrt(1.0+8.0*x))/2;
			int j = x - i*(i+1)/2;
			if (i==j){
				allHash[x] = "";
				continue;
			}
			std::string gene1 = genes[i];
			std::string gene2 = genes[j];
			int m = gene1.length(); // length of gene1
			int n = gene2.length(); // length of gene2
			int l = m+n;
			int xans[l+1], yans[l+1];
			uint64_t get_p_start = GetTimeStamp();
			penalties[probNum]=getMinimumPenalty(gene1,gene2,pxy,pgap,xans,yans);
			//printf("Time get_penalties:%1d us\n",(uint64_t)(GetTimeStamp()-get_p_start));
			// Since we have assumed the answer to be n+m long,
			// we need to remove the extra gaps in the starting
			// id represents the index from which the arrays
			// xans, yans are useful
			int id = 1;
			int a;
			for (a = l; a >= 1; a--)
			{
				if ((char)yans[a] == '_' && (char)xans[a] == '_')
				{
					id = a + 1;
					break;
				}
			}
			std::string align1="";
			std::string align2="";
			for (a = id; a <= l; a++)
			{
				align1.append(1,(char)xans[a]);
			}
			for (a = id; a <= l; a++)
			{
				align2.append(1,(char)yans[a]);
			}
			std::string align1hash = sw::sha512::calculate(align1);
			std::string align2hash = sw::sha512::calculate(align2);
			allHash[x] = sw::sha512::calculate(align1hash.append(align2hash));
			
			probNum++;		
			// Uncomment for testing purposes
/*			 std::cout << penalties[probNum] << std::endl;
			 std::cout << align1 << std::endl;
			 std::cout << align2 << std::endl;
			 std::cout << std::endl;
*/
			//#pragma omp ordered
	}
	for (x=0;x<nmax;x++){
		if (allHash[x]!=""){
			alignmentHash = sw::sha512::calculate(alignmentHash.append(allHash[x]));	
		}
	}
	return alignmentHash;
}

// function to find out the minimum penalty
// return the minimum penalty and put the aligned sequences in xans and yans
int getMinimumPenalty(std::string x, std::string y, int pxy, int pgap, int *xans, int *yans)
{
	
	int i, j; // intialising variables

	int m = x.length(); // length of gene1
	int n = y.length(); // length of gene2
	
	// table for storing optimal substructure answers
	int **dp = new2d (m+1, n+1);
	size_t size = m + 1;
	size *= n + 1;
	memset (dp[0], 0, size);

	// intialising the table
	for (i = 0; i <= m; i++)
	{
		dp[i][0] = i * pgap;
	}
	for (i = 0; i <= n; i++)
	{
		dp[0][i] = i * pgap;
	}

	// calcuting the minimum penalty
	//#pragma omp parallel for num_threads(4) shared(dp,pxy,pgap,x,y) private(i,j)
/*	
	for (i = 1; i <= m; i++)
	{
		//#pragma omp parallel for num_threads(4)
		for (j = 1; j <=n; j++)
		{	
			if (x[i - 1] == y[j - 1])
			{
				dp[i][j] = dp[i - 1][j - 1];
			}
			else
			{
				dp[i][j] = min3(dp[i - 1][j - 1] + pxy ,
						dp[i - 1][j] + pgap ,
	i					dp[i][j - 1] + pgap);
			}
		}
	}
*/	
	int diagonals,length, tile,ii, jj, iii, jjj, diff_i,diff_j,diag_i,diag_j,imax,jmax;
	int height = n;
	int width = m;
	int cores =8;// numa_num_configured_cpus();
	int di=128, dj=128;
	//int di = MAX(m/(4*cores)+1,128), dj =MAX(n/(4*cores)+1,128);
	//printf("di:%d,dj%d\n",di,dj);
	i =0, j=0;
	diagonals = (height/di) + (width/dj);
	diagonals += (height%di>0)?1:0;
	diagonals += (width%dj>0)?1:0;	
	for (int d=0; d<diagonals; ++d){
		diff_i = height-i-1;
		diff_j = j;
		diag_i = 1+((diff_i)/di);
		diag_j = 1+((diff_j)/dj);
		length = MIN(diag_i, diag_j);
		//printf("length:%d n:%d,m:%d,i:%d,j:%d\n",length,n,m,i,j);			
		#pragma omp parallel for schedule(dynamic) num_threads(MIN(length,cores)) proc_bind(spread) shared(i,j,dp,width,length,di,dj,x,y,pxy,pgap) private(iii,jjj,ii,jj,imax,jmax,tile)
		for (tile = 0; tile<length; ++tile){
			//printf("thread:%d/%d\n",omp_get_thread_num(),omp_get_num_threads());
			numa_run_on_node(0);
			//int c = sched_getcpu();
			//printf("p%d ",omp_get_place_num());
			//printf("node:%d,puid:%d\n",numa_node_of_cpu(c),c);
			ii = i + (tile*di);
			jj = j - (tile*dj);
			imax = MIN(ii+di, height);
			jmax = MIN(jj+dj, width);
			for (jjj=jj; jjj<jmax;++jjj){
				for (iii=ii;iii<imax;++iii){
					if (x[jjj] == y[iii]){
						dp[jjj+1][iii+1]=dp[jjj][iii];
					}else{
						dp[jjj+1][iii+1]=min3(dp[jjj][iii]+pxy,
								dp[jjj][iii+1]+pgap,
								dp[jjj+1][iii]+pgap);
					}
				}
			}
		}
		j+=dj;
		if (j>=width){
			j=j-dj;
			i+=di;
		}
	}
	// Reconstructing the solution
	int l = n + m; // maximum possible length
	
	i = m; j = n;
	
	int xpos = l;
	int ypos = l;
	
	while ( !(i == 0 || j == 0))
	{
		if (x[i - 1] == y[j - 1])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j - 1] + pxy == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)y[j - 1];
			i--; j--;
		}
		else if (dp[i - 1][j] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)x[i - 1];
			yans[ypos--] = (int)'_';
			i--;
		}
		else if (dp[i][j - 1] + pgap == dp[i][j])
		{
			xans[xpos--] = (int)'_';
			yans[ypos--] = (int)y[j - 1];
			j--;
		}
	}
	while (xpos > 0)
	{
		if (i > 0) xans[xpos--] = (int)x[--i];
		else xans[xpos--] = (int)'_';
	}
	while (ypos > 0)
	{
		if (j > 0) yans[ypos--] = (int)y[--j];
		else yans[ypos--] = (int)'_';
	}

	int ret = dp[m][n];

	delete[] dp[0];
	delete[] dp;
	
	return ret;
}
