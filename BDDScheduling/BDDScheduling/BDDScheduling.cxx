using namespace std;

#include<iostream>
#include<fstream>
#include<string>
#include <sys/time.h>
#include<math.h>
#include "bdd.h"
#include "bvec.h"



bool display_debug; //output the debug info.
bool display_least; //output the least info.

bool writeFile; //output the backtracking info to the outfile result.txt


//**********The following parameters should be revised correspondingly before applying to a new net!*******************
/*The h function to be used is decided by the parameters in Search() in Main(), e.g., Search("Astar",0, 0) and optimal_f=Search("Astar",3, h3_0);.*/


/*
char source_net[50]="new4x3_1111";  // the source net
#define NP 31  //number of places!
#define NT 24  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={6,11,16,23};//for h3=e_Maxi, j is started from 1. pj is the set of places which use the same uni-resource and whose total time is the biggest.
int h3_0=16; //for h3=e_Maxi  16=4+5+5+2
int P_R[]={29,30,31}; //resource places, start from 1.
*/


/*
char source_net[50]="new4x3_2111";  // the source net
#define NP 31  //number of places!
#define NT 24  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={6,11,16,23};//for h3=e_Maxi, j is started from 1. pj is the set of places which use the same uni-resource and whose total time is the biggest.
int h3_0=20; //for h3=e_Maxi  =4*2+5+5+2
int P_R[]={29,30,31}; //resource places, start from 1.
*/


/*
char source_net[50]="new4x3_2211";  // the source net
#define NP 31  //number of places!
#define NT 24  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={6,11,16,23};//for h3=e_Maxi, j is started from 1. pj is the set of places which use the same uni-resource and whose total time is the biggest.
int h3_0=25; //for h3=e_Maxi  =4*2+5*5+5+2
int P_R[]={29,30,31}; //resource places, start from 1.
*/


/*
char source_net[50]="new4x3_2221";  // the source net
#define NP 31  //number of places!
#define NT 24  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={6,11,16,23};//for h3=e_Maxi, j is started from 1. pj is the set of places which use the same uni-resource and whose total time is the biggest.
int h3_0=30; //for h3=e_Maxi  =4*2+5*5+5*5+2
int P_R[]={29,30,31}; //resource places, start from 1.
*/

/*
char source_net[50]="new4x3_2222";  // the source net
#define NP 31  //number of places!
#define NT 24  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={6,11,16,23};//for h3=e_Maxi, j is started from 1. pj is the set of places which use the same uni-resource and whose total time is the biggest.
int h3_0=32; //for h3=e_Maxi  =4*2+5*5+5*5+2*2
int P_R[]={29,30,31}; //resource places, start from 1.
*/




/*
char source_net[50]="ChenFig511";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2,5,11,13};//for h3=e_Maxi, j is started from 1.
int h3_0=16; //for h3=e_Maxi, 16=5+4+3+4
*/

/*
char source_net[50]="ChenFig522";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=32; //for h3=e_Maxi, 32=16*2
*/

/*
char source_net[50]="ChenFig533";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=48; //for h3=e_Maxi, =16*3
*/

/*
char source_net[50]="ChenFig544";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=64; //for h3=e_Maxi, =16*4
*/


/*
char source_net[50]="ChenFig555";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=80; //for h3=e_Maxi, =16*5
*/

/*
char source_net[50]="ChenFig566";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=96; //for h3=e_Maxi, 32=16*6
*/

/*
char source_net[50]="ChenFig577";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=112; //for h3=e_Maxi, 32=16*7
*/

/*
char source_net[50]="ChenFig588";  // the source net
#define NP 21  //number of places!
#define NT 14  //number of transitions!
#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 8 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
const bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
int h3pj[]={2, 5, 11, 13};//for h3=e_Maxi, j is started from 1.
int h3_0=128; //for h3=e_Maxi, 32=16*8
*/





char source_net[50]="Chen2011Big11112";  // the source net
bool isIncidenceMatrix = false;  //if not transposed, then true; if transposed, then false.
#define NP 53  //number of places!
#define NT 38  //number of transitions!
#define MAXNUMBER 1023 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
#define MAXBITS 10 // the number of bits representing the above information. MAXBITS=RoundUp(log_2(MAXNUMBER+1))
int h3pj[]={2, 17, 34, 36};//for h3=e_Maxi, j is started from 1.!!!!!!!!!!!!!
int h3_0=32; //for h3=e_Maxi!!!!!!!!!!!
int P_R[]={37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48}; //resource places, start from 1.



//*********************************************************************************************************************

const int number_var_bigM = 3*NP+1;




bdd open;
bdd closed;

bdd et[NT];
bdd haveRemainTime[NT];

char resultFile[50];
ofstream outfile;
ifstream infile;

int initial_marking[NP];  //the initial marking
int final_marking[NP]; //the goal marking
int D[NP];  //processing time for all places

int rM[NP][NT]; //incidence matrix
int rIM[NT][NP]; //transposed incidence matrix (The code uses rIM, not rM!)

int loopNo;




bvec F;
bvec M[NP];  //M[i] = bvec_varfdd(i); //P82; after using fdd_extdomain();  //bvec* M[NP]??
bvec R1[NP];
bvec R2[NP]; // It can handle the nets whose operation places have at most two tokens.



struct schedulingState
{
	int f;
	int mu[NP]; //represents \mu in the pseudo code of the paper
	int r1[NP]; //the remaining time of the fisrt token in each place
	int r2[NP]; //the remaining time of the second token in each place
}; //do not forget the semicolon



struct stateD
{
	schedulingState s; //scheduling state
	int d; //decreased time to reversely reach the state
}; //do not forget the semicolon

stateD sd; //BddToAStateHandler() will change sd.


stateD BackT(bdd, schedulingState, schedulingState, int);
bdd Constraints(schedulingState, int, int, int);
bool IsInitial(schedulingState);
bool IsResource(int); //decide whether p_i+1 is a resource place





void InitMarking()  //Constructs initial_marking, D, and final_marking from the input outfile xxx_init.txt
{
	char inputFile[50];
	strcpy(inputFile, source_net);
	strcat(inputFile, "_init.txt");			
	ifstream fin(inputFile);
        string str_init;
        int i=0;
        while(fin >> str_init && i<3*NP)
        {
		if(i<NP)
		        initial_marking[i]=atoi(str_init.c_str());		 	
		else if(i>=NP && i<2*NP)		
			D[i-NP]=atoi(str_init.c_str());				
		else
			final_marking[i-2*NP]=atoi(str_init.c_str());	
		i++;
        }
}


void InitMatrix() //When the input outfile contains a incidence matrix, constructs incidence matrix rM from the input outfile xxx_matrix.txt
{
	char inputFile[50];
	strcpy(inputFile, source_net);
	strcat(inputFile, "_matrix.txt");
	ifstream fin(inputFile);
        string str_data;
        int i=0,j=0;
        while(fin >> str_data)
        {
                if(i == NP)
                {
                        break ;
                }
                rM[i][(j-i*NT)]=atoi(str_data.c_str());
                if(((j+1)%NT)==0 && j!=0)
                {
                        i++;
                }
                ++j;
        }
}

void TransposeMatrix()  //Constructs transposed incidence matrix rIM from rM since this program uses the transposed one rIM.
{
	for(int i=0;i<NP;i++)
	{
		for(int j=0;j<NT;j++)
		{
			rIM[j][i] = rM[i][j];
		}
	}	
}


void InitTransposedMatrix() //When the input outfile contains a transposed incidence matrix derived from INA, constructs transposed incidence matrix rIM from it.
{
	char inputFile[50];
	strcpy(inputFile, source_net);
	strcat(inputFile, "_matrix.txt");
	ifstream fin(inputFile);
        string str_data;
        int i=0,j=0;
        while(fin >> str_data)
        {
                if(j == NT)
                {
                        break ;
                }
                rIM[j][(i-j*NP)]=atoi(str_data.c_str());
                if(((i+1)%NP)==0 && i!=0)
                {
                        j++;
                }
                i++;
        }
}



bdd SubtractOne(bdd from,bdd B,int n,int a[])  //a recursive function for bit borrow.  a[] stores the positions of the bits of a place. B=bddtrue denotes subtracting 1; bddfalse not. n=0  //Note that the addition and subtraction in bits are more efficient than those in bvec. (11111)_2 and (31)_10
//At the beginning, n=0.
{
        bdd var,var1,var2,tmp1,tmp2,f1,f2,b1,b2,r;

        if(from == bddfalse)
                return bddfalse;        
        if(B == bddfalse)  //if need not to substract 1.
                return from;        

        if(n == MAXBITS-1) //if the calculation is being applied to the highest bit
        { 
                var = bdd_ithvar(a[n]);
                tmp1 = bdd_constrain(from,var);
                f1 = bdd_and(bdd_not(var),tmp1);
                return f1;
        }
        else
        {
                var = bdd_ithvar(a[n]);
                tmp1 = bdd_constrain(from,var);		
                if(tmp1 != bddfalse)
                {
                        f1 = bdd_and(bdd_not(var),tmp1);
                        b1 = bddfalse;
                }
                else
                        f1 = bddfalse;
                
                tmp2 = bdd_constrain(from,bdd_not(var));
			
                if(tmp2 != bddfalse)
                {
                        f2 = bdd_and(var,tmp2);
                        b2 = bddtrue;
                }
                else
                        f2 = bddfalse;
                        
                n++;
                var1 = SubtractOne(f1,b1,n,a);
                var2 = SubtractOne(f2,b2,n,a);
                r = bdd_or(var1,var2);
                
                return r;
        }
}

bdd AddOne(bdd from,bdd B,int n,int a[])  //a recursive function for bit carry. a[] stores the positions of the bits of a place. B=bddtrue denotes adding 1; bddfalse not. n=0
{
        bdd var,var1,var2,tmp1,tmp2,f1,f2,b1,b2,r;
        if(from == bddfalse)  //if addend has nothing.
                return bddfalse;
        if(B == bddfalse)
                return from;
        
        if(n == MAXBITS-1) //if the calculation is being applied to the highest bit
        { 
                var = bdd_ithvar(a[n]);
                tmp1 = bdd_constrain(from,bdd_not(var));
                f1 = bdd_and(var,tmp1);
                return f1;
        }
        else
        {
                var = bdd_ithvar(a[n]);
                tmp1 = bdd_constrain(from,bdd_not(var));
                if(tmp1 != bddfalse)
                {
                        f1 = bdd_and(var,tmp1);
                        b1 = bddfalse;
                }
                else
                        f1 = bddfalse;
                
                tmp2 = bdd_constrain(from,var);
                if(tmp2 != bddfalse)
                {
                        f2 = bdd_and(bdd_not(var),tmp2);
                        b2 = bddtrue;
                }
                else
                        f2 = bddfalse;
                
                n++;
                var1 = AddOne(f1,b1,n,a);
                var2 = AddOne(f2,b2,n,a);
                r=bdd_or(var1,var2);

                return r;
        }
}


/*
bdd addBDD(bdd addend,bdd adding, bdd B, int a[])
{
	bdd var,var1,var2,tmp1,tmp2,f1,f2,b1,b2,r;

        if(addend == bddfalse) // if addend has nothing.
                return adding;
                
                
        if(B == bddfalse)
                return from;
                
        
}
*/



void  AllsatPrintHandler(char* varset, int size)  //size is the number of bits in each marking. Use bdd_allsat(bdd r, AllsatPrintHandler) to print all markings.
//cout<<fddset<<from<<endl and  fdd_printset(bdd) have a similar result to that of bdd_allsat(from, AllsatPrintHandler).
{
	int f_dec=0;	
	int M_dec[NP], R1_dec[NP], R2_dec[NP];
	for(int i=0;i<MAXBITS;i++)
	{
		f_dec += varset[i*number_var_bigM]*(int)pow(2,i);		
	}
	for(int n=0;n<NP;n++)
	{
		M_dec[n]=R1_dec[n]=R2_dec[n]=0;
		for(int i=0;i<MAXBITS;i++)
		{
			M_dec[n] += varset[n+1+i*number_var_bigM]*(int)pow(2,i);
			R1_dec[n] += varset[n+NP+1+i*number_var_bigM]*(int)pow(2,i);
			R2_dec[n] += varset[n+2*NP+1+i*number_var_bigM]*(int)pow(2,i);
		}
	}

	outfile<<"<f="<<abs(f_dec)<<", M=(";	
	for(int n=0;n<NP-1;n++)		
		outfile<<abs(M_dec[n])<<", ";
		
	outfile<<abs(M_dec[NP-1])<<"), R1=(";
	for(int n=0;n<NP-1;n++)		
		outfile<<abs(R1_dec[n])<<", ";
		
	outfile<<abs(R1_dec[NP-1])<<"), R2=(";
	for(int n=0;n<NP-1;n++)		
		outfile<<abs(R2_dec[n])<<", ";
				
	outfile<<abs(R2_dec[NP-1])<<")>"<<endl;	
	
}




void  BddToAStateHandler(char* varset, int size)  //size is the number of bits in each marking. Use bdd_allsat(bdd r, BddToAStateHandler) to transform a bdd r to sd.s. Note that use r=bdd_fullsatone(found) first to ensure there is only one state in the bdd.
{
	int f_dec=0;	
	int M_dec[NP], R1_dec[NP], R2_dec[NP];
	for(int i=0;i<MAXBITS;i++)
	{
		f_dec += varset[i*number_var_bigM]*(int)pow(2,i);		
	}
	for(int n=0;n<NP;n++)
	{
		M_dec[n]=R1_dec[n]=R2_dec[n]=0;
		for(int i=0;i<MAXBITS;i++)
		{
			M_dec[n] += varset[n+1+i*number_var_bigM]*(int)pow(2,i);
			R1_dec[n] += varset[n+NP+1+i*number_var_bigM]*(int)pow(2,i);
			R2_dec[n] += varset[n+2*NP+1+i*number_var_bigM]*(int)pow(2,i);
		}
	}
	
	sd.s.f=abs(f_dec);
	
	for(int i=0;i<NP;i++)
	{
		sd.s.mu[i]=abs(M_dec[i]);
		sd.s.r1[i]=abs(R1_dec[i]);
		sd.s.r2[i]=abs(R2_dec[i]);	
	}		
}


//for h3=e_Maxi
//decide if Place x is in the vector h3pj
bool IsInH3pj(int x)
{
	int n_h3pj=sizeof(h3pj)/sizeof(h3pj[0]); //IMPORTANT!!!!
	for(int n=0; n<n_h3pj;n++)
		if(x==h3pj[n])
			return true;
	return false;
}



bdd Img(bdd from, int hmethod) //To generate a set of submarkings of from
{
        bdd generatedMarkings = bddfalse;
	bdd markingsEnableTi; //the set of markings which enable t_i.
	bdd result_tmp;
	
	
	// et[i]: The constraint representing that t_i is enalbed. It should be combined with a marking set.
	// haveRemainTime[i]: The constraint representing there exists any remaining process time in any pre-place of t_i. It should be combined with a marking set.
	for(int i=0;i<NT;i++) //t_i     
	{
		et[i] = bddtrue;  //It is a constraint. It should be combined with a marking set.
		haveRemainTime[i] = bddfalse; //It is a constraint. It should be combined with a marking set.
		for(int j=0;j<NP;j++) //p_j
			if(rIM[i][j] < 0) //p_j->t_i
			{
				et[i] &= bvec_gte(M[j],bvec_con(MAXBITS,1)); // Note that: here, the value of an arc is always 1. For ordinary nets.
				if(D[j]!=0)  //only consider operational input places of t_i.
				{
					haveRemainTime[i] |= bvec_equ(M[j], bvec_con(MAXBITS, 1)) & (bvec_gth(R1[j],bvec_con(MAXBITS,0)) | bvec_gth(R2[j],bvec_con(MAXBITS,0))); 
					haveRemainTime[i] |= bvec_equ(M[j], bvec_con(MAXBITS, 2)) & bvec_gth(R1[j],bvec_con(MAXBITS,0)) & bvec_gth(R2[j],bvec_con(MAXBITS,0));
				}
			}
	}
        
        
        for(int i=0;i<NT;i++)  //for each transition t_i
        {   	
        	
		markingsEnableTi = from & et[i]; //et[i]: the condition under which t_i is enabled in a marking. Here, only handle the markings which are in From and can enable t_i.
				
                if(markingsEnableTi!=bddfalse)  // Here, cannot use markingsEnableTi==bddtrue
                {   

		        bdd tmp1, tmp2;
        		int a[MAXBITS];        		
        		
                	             	
                	//11111111111111111. Subtract all remaining process time in each t_i's pre-places, and add the corresponding time to "f" (f=g+h).
			//Delta is the maximum of the remaining process times of all places p\in ^{\bullet} t_i.                	
                	//Subtract delta from Mr and add it to f.  //The method is to iteratively subtract one from all Mr if any pre-post of the transition has remaining process time.
                	//Add delta to f       
                	//Handle h.
                	// (F, M, R1, R2)            	                	

		        tmp1 = markingsEnableTi & haveRemainTime[i];  //t_i
		        
		        result_tmp = markingsEnableTi - tmp1;
		        
		        
		        
		         //For all markings which have remain time to fire t_i, subtract 1 from each place whose remain process time is greater than 0 until t_i is ready to fire and add 1 to f at each iteration.
		        while(tmp1 != bddfalse) // While there exists some pre-place of t_j whose Mr > 0 (consider different tokens in a place), then subtractOneFromAllMr and addOneToF. tmp1 contains markings which have remain time in some place to fire t_i.
		        {        	    
		        	for(int j=0;j<NP;j++)  //R1 and R2 of all oprational places will be subtracted one if they are greater than 0.
		        		if(D[j]>0)
			        	{
			        					        		
			        		//subtract 1 from R1 of all operational places if their R1>0.//
			        		//Uses SubtractOne() when the minuend is greater than 0, Mr cannot be negative.
			        		bdd tmp1_gt0 = tmp1 & bvec_gth(R1[j],bvec_con(MAXBITS,0));
			        		
			        		
			        		if (tmp1_gt0!=bddfalse) 
			        		{
			        			tmp1 = tmp1-tmp1_gt0;
			        			
				        		for(int s=0;s<MAXBITS;s++)        
		                                		a[s]=j+1+NP+s*number_var_bigM;       //R1[j]   
		                                	tmp1_gt0 = SubtractOne(tmp1_gt0,bddtrue,0,a);//R1[j]-- 
		                                	
		                                	if(hmethod==3 && IsInH3pj(j+1))//f(M')=f(M)+c(M,M')+h(M')-h(M),here perform h(M')-h(M) for h3=e_Maxi
		                                	{
		                                		// Subtract 1 from f.		
						        	for(int s=0;s<MAXBITS;s++)        
				                                	a[s]= s*number_var_bigM;  //f                        
				                                tmp1_gt0 = SubtractOne(tmp1_gt0,bddtrue,0,a);  //f--
		                                	}
		                                	
		                                	tmp1 |= tmp1_gt0;
	                                	}
	                                	

	                                	
	                                	tmp1_gt0 = tmp1 & bvec_gth(R2[j],bvec_con(MAXBITS,0));
			        		
			        		
	                                	if (tmp1_gt0 != bddfalse)
			        		{
			        			tmp1 = tmp1-tmp1_gt0;
			        			
		                                	for(int s=0;s<MAXBITS;s++)        
		                                		a[s]=j+1+2*NP+s*number_var_bigM;        //R2[j]                        
		                                	tmp1_gt0 = SubtractOne(tmp1_gt0,bddtrue,0,a);   //R2[j]--
		                                	
		                                	if(hmethod ==3 && IsInH3pj(j+1))//f(M')=f(M)+w(M,M')+h(M')-h(M),here perform h(M')-h(M) for h3=e_Maxi
		                                	{
		                                		// Subtract 1 from f.
						        	for(int s=0;s<MAXBITS;s++)        
				                                	a[s]= s*number_var_bigM;  //f                        
				                                tmp1_gt0 = SubtractOne(tmp1_gt0,bddtrue,0,a);  //f--
		                                	}
		                                	
		                                	tmp1 |= tmp1_gt0;
		                                }      
		                                
			        	}// for(int j=0;j<NP;j++)          	
		        	
		        	
                                      	
		        	// Add 1 to f.
		        	for(int s=0;s<MAXBITS;s++)        
                                	a[s]= s*number_var_bigM;  //f                        
                                tmp1 = AddOne(tmp1,bddtrue,0,a);  //f++
                                
					
                                tmp2 = tmp1 & haveRemainTime[i];  //t_i
                                result_tmp |= (tmp1 - tmp2);  
                                tmp1 = tmp2;         
                                
                                		        	
		        }//while(tmp1 != bddfalse)     
		               
		       
			if (hmethod==1) //h=t_depth-1; h_0=12
	                	//...........f(M)+h(M')-h(M); h=depth-1 when the fired t is a post-t of an operational place. 
	                	if (i==2 ||i==4 ||i==6 ||i==8 ||i==10 ||i==12 ||i==14 ||i==16 ||i==18 || i==20 ||i==22 ||i==24)
	                	{
	                		int a[MAXBITS];  //the value is the number of bits representing a place
	                                for(int s=0;s<MAXBITS;s++)        
	                                        	a[s]=s*number_var_bigM;                                   
	                                                
	                                result_tmp = SubtractOne(result_tmp,bddtrue,0,a);
	                	}
                	
                	
                	
                	// 2222222222222222. Move tokens in M. 
         		for(int j=0;j<NP;j++)  //p_j
         		{ 
                                if(rIM[i][j] < 0)   //rIM: transposed incidence matrix; if some tokens should be taken away from some input places by the firing of the transition. Note that it only handles the arcs with value 1.
                                {  
                                        int a[MAXBITS];  //the value is the number of bits representing a place
                                        for(int s=0;s<MAXBITS;s++)        
                                        	a[s]=j+1+s*number_var_bigM; //2017/06/17 - 09:58 - Bo Huang  The positions of bits of tokens in a place.                                 

						//Note that: fdd_extdomain(int* dom, int num) extends the set of finite domain blocks with the NP domains in dom. Each domain allocates log_2(dom[i]) BDD variables. The ordering is interleaved.
                                        result_tmp = SubtractOne(result_tmp,bddtrue,0,a);//result_tmp is the operated BDD; bddtrue denotes to subtract 1; "a" stores the positions of bits of places. Handle all markings in a set. 
                                        
                                        //2017/06/14 - 19:21 - Bo Huang   Uses bvec_add() and bvec_sub(). They can automatically carry and borrow bits.
                                        // However, operations on bits are more efficient than on bvec.
                                            
                                }
                                if(rIM[i][j] > 0) // if some tokens should be added to the place via the arc
                                { 
                                        int a[MAXBITS];
                                        for(int s=0;s<MAXBITS;s++)
                                        	a[s]=j+1+s*number_var_bigM;
               
					result_tmp = AddOne(result_tmp,bddtrue,0,a);					
                                }
                        }
                        
                	
                	//333333333333333: Add the process times to the R1 or R2 of all output places of t_i. Should consider the number of tokens in a place.
                	//The current version is for the number of tokens in any operation place is at most 2. It can be revised for other situations.
                	bdd R,tmp3;
                	bdd result_3=bddfalse;
                	
                	
                	for(int j=0;j<NP;j++) //p_j
                	{         
		         	if(rIM[i][j] > 0 && D[j]!=0) //t_i-> p_j and p_j is an operation place.
		         	{
		         		
		         		//Case 1: R1[j]<=0.
		         		R = bvec_lte(R1[j], bvec_con(MAXBITS, 0));  //constraint: R1[j]<=0
		         		tmp1 = result_tmp & R;		 //constraint: R1[j]<=0
		         		
		         		//Case 2: R1[j]>0.
		         		R = bvec_gth(R1[j], bvec_con(MAXBITS, 0));  //constraint: R1[j]>0
		         		tmp2 = result_tmp & R;		//constraint: R1[j]>0
		         		
		         		
		         		if(tmp1 != bddfalse)
		         		{
		         			result_tmp = result_tmp - tmp1;
		         			
		         			tmp1=bdd_exist(tmp1, fdd_ithset(j+1+NP)); //P34: Remove all occurences in 1# of the variables in 2# by existential quantification.
		         			tmp1 &= fdd_ithvar(j+1+NP,D[j]);//Then R1[j]=D[j]
		         			
		         			result_tmp |= tmp1;		         			
		         		}  
		         		
		         		
		         		if(tmp2 != bddfalse)
		         		{
		         			result_tmp = result_tmp - tmp2;
		         			
		         			tmp2=bdd_exist(tmp2, fdd_ithset(j+1+2*NP));
		         			tmp2 &= fdd_ithvar(j+1+2*NP,D[j]); //Then R2[j]=D[j]
		         			
		         			result_tmp |= tmp2;		         			
		         		} 
		         		
		         	}//if(rIM[i][j] > 0 && D[j]!=0)
		         }//for(int j=0;j<NP;j++) //p_j	 
		             
                	                	
                        generatedMarkings |= result_tmp;
                        
                        
                }//if(markingsEnableTi!=bddfalse) 
                
        }//for(int i=0;i<NT;i++)   t_i
        
        return generatedMarkings;
}




int Traverse(string choice, int hmethod)  // Search in the state space.  If it finds a path with Astar algorithm, then return the optimal f value; if it finds a path with BFS, return 0; if no result, return -1.
{
	
	bdd R_Mt =bddtrue;  //M in bigM is restricted to the final marking.	
	for(int i=0;i<NP;i++)
		R_Mt  &= fdd_ithvar(i+1,final_marking[i]);  //P86, fdd_ithvar() returns the BDD that defines the value 2# for the finite domain block 1#.	
	
	loopNo=1;  //Global variable for debugging.	
	
	bdd from, R, found;
	bdd succ = bddfalse ; 	
	
	
	//Breadth First Search
	if(choice=="BFS" || choice=="bfs")
		while(1)
		{
			if(display_least==false)
			{
				cout<<"=====================loopNo:"<<loopNo<<"====================="<<endl;
				cout<<"Number of markings in open: "<<bdd_satcount(open)<<endl;		
				if(display_debug==true)
				{
					bdd_allsat(open, AllsatPrintHandler);
					cout<<endl;					
				}
			}
			
			if(open == bddfalse)
		        {
		                cout<<"Open is NULL! No result for this search!!"<<endl;
		                return -1;               
		        } 	
	
			found = open & R_Mt; //Note that: should not use succ & R_Mt!!
			if(found!= bddfalse)  //Finds the markings whose marking is Mt, i.e., if open contains Mt; open consists of 2N*2 interleaved FDD variables.
			{
			                	cout<<"A result has been found for the NET: "<<source_net<<"."<<endl;
			                	cout<<"Format of fdd_printset(found): <0:f, 1-"<<NP<<":M, "<<NP+1<<"-"<<2*NP<<":R1, "<<2*NP+1<<"-"<<3*NP<<":R2>"<<endl;
			                	cout<<fddset<<found<<endl;
			                	
			                	
			                	if(writeFile==true)
			                	{
				                	outfile<<"A result has been found for the NET: "<<source_net<<"."<<endl;
				                	outfile<<"Format of fdd_printset(found): <0:f, 1-"<<NP<<":M, "<<NP+1<<"-"<<2*NP<<":R1, "<<2*NP+1<<"-"<<3*NP<<":R2>"<<endl;
				                	outfile<<fddset<<found<<endl;
			                	}

			                	
			                	
			                	
			                	if(display_least==false)
			                	{		               		
			                		cout<<"Mt(bdd_allsat)=";
			                		bdd_allsat(found, AllsatPrintHandler); 
			                		cout<<endl;
			                	
			                	
			                		cout<<"----------open------bdd_satcount(open):"<<bdd_satcount(open)<<endl;
			                		fdd_printset(open);
			                		cout<<endl;
			                		cout<<"----------found-----bdd_satcount(found):"<<bdd_satcount(found)<<endl;
			                		fdd_printset(found);
			                		cout<<endl;
			                	}

			                	
			                	
			                    	return 0;               
			} 
			
			
			closed |= open; //To compute the number of explored states.
			open = Img(open, hmethod);  // At the beginning, from only contains M0.
			
			loopNo++;         

		}//while(1)	
	
	else if(choice=="Astar" || choice=="astar")
	{
		int least_f=0; //since we use a consistent heuristic, so f is non-decreasing.
		while(1)
		{
			if(open == bddfalse)
		        {
		                cout<<"Open is NULL! No result for this search!!"<<endl;
		                return -1;               
		        } 		        		        
		
			for(int f_scan=least_f;;f_scan++)
			{					        
			        if(display_least==false)
			        {
				        cout<<"=====================loopNo:"<<loopNo<<", f_scan: "<<f_scan<<"====================="<<endl;
					cout<<"Number of markings in OPEN: "<<bdd_satcount(open)<<endl;		
					cout<<"Number of markings in CLOSED: "<<bdd_satcount(closed)<<endl;	
				}				
				
				R = fdd_ithvar(0, f_scan);		
				from = open & R;				
				
				if(from != bddfalse)
				{
					least_f=f_scan;  //since we use a consistent heuristic, so f is non-decreasing.
					
					found = from & R_Mt; //Note that: we should not use succ & R_Mt!!
			        	if(found!= bddfalse)  //Finds the states whose marking is Mt, i.e., if open contains Mt; open consists of 2N*2 interleaved FDD variables.
			        	{			
				     		
			                	
			                	cout<<"A result has been found for the NET: "<<source_net<<"."<<endl;
			                	cout<<"Format of fdd_printset(found): <0:f, 1-"<<NP<<":M, "<<NP+1<<"-"<<2*NP<<":R1, "<<2*NP+1<<"-"<<3*NP<<":R2>"<<endl;
			                	cout<<fddset<<found<<endl;
			                	
			                	if(writeFile==true)
			                	{
				                	outfile<<"A result has been found for the NET: "<<source_net<<"."<<endl;
				                	outfile<<"Format of fdd_printset(found): <0:f, 1-"<<NP<<":M, "<<NP+1<<"-"<<2*NP<<":R1, "<<2*NP+1<<"-"<<3*NP<<":R2>"<<endl;
				                	outfile<<fddset<<found<<endl;
			                	}
			                	
			                	
			                	if(display_least==false)
			                	{
			                		cout<<"Mt(bdd_allsat)=";
			                		bdd_allsat(found, AllsatPrintHandler);
			                		cout<<endl;
			                	}		                
			                	
			                    	return least_f;               
			         	} 		
					
					
					
					if(display_debug==true && display_least==false)
					{
						cout<<"loopNo:"<<loopNo<<"--------"<<endl;
						cout<<"Number of markings in from: "<<bdd_satcount(from)<<endl;						
						bdd_allsat(from, AllsatPrintHandler);
					}
					
					
					succ = Img(from, hmethod);  // At the beginning, from only contains M0.

					
					
					closed |= from;	
					succ=succ-closed; 
					
					if(display_debug==true && display_least==false)
					{
						cout<<"loopNo:"<<loopNo<<"--------"<<endl;
						cout<<"Number of new markings in succ: "<<bdd_satcount(succ)<<endl;
						bdd_allsat(succ, AllsatPrintHandler);
						cout<<endl;
					}
					open = open - from;
					open |= succ;
					
					from = open & R; //Use the same f_scan. if it generates new markings with the same f value.
					
					loopNo++;
					
					break;  //Break the FOR structure. Only if from==bddfalse, then f_scan++ in FOR.				     
			        }//if(from != bddfalse)	    			        			     

			}//for(int f_scan=least_f;;f_scan++)
		}//while(1)
	}//else if(choice=="Astar" || choice=="astar")
	else
	{
		cout<<"The parameter of hmethod is wrong!"<<endl;
		return -1;
	}
}




int Search(string choice, int hmethod, int h0)  // Obtians an optimal path from the initial marking to the given final marking via A* search
{
	int i ;
	
	static int dom[number_var_bigM]; //const int number_var_bigM = 3*NP+1; (F,M,R1,R2)
	for(int i=0;i<number_var_bigM;i++)
		dom[i]= MAXNUMBER; //#define MAXNUMBER 255 // the maximal number for the possible makespan, tokens in a place, remaining time in a place, and the number of a transition.
	
	
	bdd_init(10000000,1000000);
	fdd_extdomain(dom,number_var_bigM);  //P84, Extends the set of finite domain blocks with the NP domains in dom. Each domain allocates log_2(dom[i]) BDD variables. The ordering is interleaved.
	
		
	
	// renames the elements in bigM=(F,M,R1,R2).  
	F=bvec_varfdd(0);
	for(i=0;i<NP;i++)
		M[i] = bvec_varfdd(i+1); //bvec M[NP]  // P82,bulds a Boolean vector which will include exactly the variables used to define the FDD variable block i.
	for(i=0;i<NP;i++)
		R1[i] = bvec_varfdd(i+NP+1);
	for(i=0;i<NP;i++)
		R2[i] = bvec_varfdd(i+2*NP+1);
		
	
	bdd M0_big=bddtrue; // the initial marking (F,M,R1,R2)
	M0_big &= fdd_ithvar(0,h0);  // f(M_0)=h0;
		
	for(i=0;i<NP;i++) // M(M_0)=the initial marking; R1(M_0)=0; R2(M_0)=0;
	{
		M0_big &= fdd_ithvar(i+1,initial_marking[i]); //P86  Returns the BDD that defines the value initial_marking[i] for the finite domain block i+1. Do not operate a bdd with a fdd by using &=.
		M0_big &= fdd_ithvar(i+NP+1,0); //R1
		M0_big &= fdd_ithvar(i+2*NP+1,0); //R2
	}
		
	
	
	cout<<"cout<<fddset<<M0_big:"<<endl;
	cout<<fddset<<M0_big<<endl;
	cout<<endl;
	
	if(writeFile==true)
	{
		outfile<<"M0_big (bdd_allsat):"<<endl;
		bdd_allsat(M0_big, AllsatPrintHandler);
		outfile<<endl;
	}

	
		
	open = bddfalse;
	open |= M0_big;  // At the beginning, only M0 is in open.
	
	if(writeFile==true)
	{
		outfile<<"Number of scheduling states in OPEN at the beginning: "<<bdd_satcount(open)<<endl;	
		outfile<<"Number of scheduling states in OPEN at the beginning (logarithm in base 2): "<<bdd_satcountln(open)<<endl;
	}

	return Traverse(choice, hmethod);
	
}

int ite_No; //ite_No represents e of the pseudo code.

//Looks for an optimal path from s_G to s_0 in CLOSED and print it.
void Construct(bdd closed, schedulingState s_0, schedulingState s_G)
{
	int path_length=100;
	int path[path_length][2];
	
	//stateT St_pre;
	//initialization
	for(int i=0;i<path_length;i++) 
	{
		for(int j=0;j<2;j++)
		{
			path[i][j]=0;
		}
	}	
	
	ite_No=0; //ite_No represents e in the pseudo code.
	schedulingState s_p=s_G; //s_p represents s' in the pseudo code.
	int g_reversed=0;
	
	
	while(IsInitial(s_p)==false) //if it does not reach the initial state.
	{
		cout<<"Iteration: "<<ite_No<<endl;
		if(writeFile==true)
		{		
			outfile<<endl;
			outfile<<"=========================================================================="<<endl;
			outfile<<"Iteration: "<<ite_No<<endl;
			outfile<<"Current state (s): ";
			outfile<<"f="<<s_p.f<<endl;
			
			for(int i=0;i<NP;i++)
			{
				if(s_p.mu[i]>0)
					outfile<<"mu("<<i+1<<"):"<<s_p.mu[i]<<", ";
			}
			for(int i=0;i<NP;i++)
			{
				if(s_p.r1[i]>0)
					outfile<<"r1("<<i+1<<"):"<<s_p.r1[i]<<", ";
			}
			for(int i=0;i<NP;i++)
			{
				if(s_p.r2[i]>0)
					outfile<<"r2("<<i+1<<"):"<<s_p.r2[i]<<", ";
			}
			outfile<<endl;
		}
	
	
		schedulingState s=s_p;
		
		//Looks for the transitions revsersly enabled at s (considering Ri) 
		bool rEt[NT]; //whether t is reversely enabled at marking s_G
		
		if(writeFile == true)        
			outfile<<"Reversely enabled transitions at the current state (s): ";
			
		for(int j=0; j<NT; j++)
		{
			rEt[j]=true;	
			for(int i=0; i<NP; i++)
			{
				if(rIM[j][i]>0 && s.mu[i]<rIM[j][i] || rIM[j][i]>0 && s.mu[i]>=rIM[j][i] && s.r1[i]<D[i] && s.r2[i]<D[i]) //the arc weight is 1
				{
					rEt[j]=false;
					break;
				}
				
			}
			if(rEt[j]==true && writeFile == true)
				outfile<<j+1<<"  ";			
						
		}
		if(writeFile == true) 		
			outfile<<endl<<endl;	
		
		for(int j=0;j<NT;j++)
		{
			if(rEt[j]==true)
			{
				stateD sd=BackT(closed,s_0,s,j); //the index of j begins with 0; 
				s_p=sd.s;
				if(s_p.f!=-1)  //if a pre-state is found in Closed.
				{			
					path[ite_No][0]=s_G.f-g_reversed;
					path[ite_No][1]=j;	
					if(writeFile == true) 
						outfile<<"path["<<ite_No<<"]: g="<<path[ite_No][0]<<", t="<<path[ite_No][1]+1<<endl;
					g_reversed+=sd.d;
					ite_No++;
					break; 
				}

			}
		}
	
	
	}//while(IsInitial(s_p)==false) 
	
	
	
	//print the path
	cout<<endl;
	cout<<"The scheduling path:"<<endl;
	for(int i=ite_No-1;i>=0;i--)
		cout<<ite_No-i<<": (g="<<path[i][0]<<", t=t_"<<path[i][1]+1<<");"<<endl;
	
	
	if(writeFile==true)
	{
		outfile<<endl;
		outfile<<"The scheduling path:"<<endl;
		for(int i=ite_No-1;i>=0;i--)	
			outfile<<ite_No-i<<": (g="<<path[i][0]<<", t=t_"<<path[i][1]+1<<");"<<endl;			

	}
	
	return;
	
}

//Decide whether s equals the initial state
bool IsInitial(schedulingState s)
{
	for(int i=0;i<NP;i++)
		if(s.mu[i]!=initial_marking[i])
			return false;
	return true;
}


//Decide whether p_i+1 is a resource place
bool IsResource(int i)
{
	int count_resource=sizeof(P_R)/sizeof(P_R[0]); //IMPORTANT!!!!
	for(int k=0; k<count_resource;k++)
		if(i+1==P_R[k]) //P_R[i] start from 1.
			return true;
	return false;
}


// Reversely firing t from (s, t_{t+1}) to s' in CLOSED. d is the decreased time to reach s'. Note that s' depends not only t but also t' which is reversely enabled at s'. s'.f=-1 means no satisfied state in Closed.
//Note that s'.mu and s'.f are decided by t, but s'.r is also decided by t' which is a transition reversely enabled at s'!!!
stateD BackT(bdd closed, schedulingState s_0, schedulingState s, int t)  // t starts from 0.
{

	
	//Token moving for the reversly firing of t in all places
	for(int i=0;i<NP;i++)
	{
		if(rIM[t][i]>0)
			s.mu[i]=s.mu[i]-rIM[t][i]; //use the transposed incidence matrix rIM
		if(rIM[t][i]<0)
			s.mu[i]=s.mu[i]-rIM[t][i];
	}
	
	
	if(IsInitial(s)==true)
	{
		sd.s=s_0;
		sd.d=0;  
		return sd;
	}
	
	
	
	//Handle the remaining time of tokens in the output activity place of t. At most one such place.
	for(int i=0;i<NP;i++)
	{
		if(rIM[t][i]>0 && D[i]>0)
		{
			if(s.r2[i]==D[i])
				s.r2[i]=0;
			else
				s.r1[i]=0;
		}
	}
	
	//Handle r1 and r2 of other activity places
	schedulingState s_p=s; //s_p denotes s'

	
	
	bool rEt[NT]; //whether t is reversely enabled at marking s_p
	
	
	if(writeFile==true)
	{
		outfile<<endl;
		outfile<<"-------------------------------------------------------------------------"<<endl;
		outfile<<"The selected enabled transition at the current state (s):"<<endl;
		outfile<<"t="<<t+1<<endl;
		outfile<<"Reversely enabled transitions at the marking of preceding state (s'.mu): "<<endl;
		outfile<<"t'=";	
		
	}
		
	for(int j=0; j<NT; j++)
	{
		rEt[j]=true;	
		for(int i=0; i<NP; i++)
			if(rIM[j][i]>0 && s_p.mu[i]<rIM[j][i])
			{
				rEt[j]=false;		
				continue;
			}	
		if(rEt[j]==true && writeFile==true)
			outfile<<j+1<<"  ";				
	}
	if(writeFile==true)
	{
		outfile<<endl;
		outfile<<"-------------------------------------------------------------------------"<<endl;
	}
	
	for(int j=0;j<NT;j++)
	{		
		if(rEt[j]==true)
		{
			bdd C=bddtrue; //note that: not to use true or bddfalse!!!
			
			//Constraints of M.
			for(int i=0;i<NP;i++)
				C &= bvec_equ(M[i], bvec_con(MAXBITS,s_p.mu[i]));
			
			
			int d,p_p,p_pp; //d is the time required for p_p to reversely enable t_{j+1} at s_p; p_p denotes p' and p_pp denotes p'' in the psuedo code of the paper.
			if(writeFile==true)
				outfile<<"t'=t_"<<j+1<<endl;
				
			for(int i=0;i<NP;i++)
			{
				if(rIM[j][i]>0 && initial_marking[i]==0) //initial_marking[i]=0	 means P_A+end places. In the reversely firing, end places can be excluded.  D[i] are allowed to be 0. 
				{
					p_p=i;
					break; //Only one such place p_{i+1}.
				}
			}
			if(writeFile==true)			
				outfile<<"p'=p_"<<p_p+1<<endl;
			
			
			for(int i=0;i<NP;i++)
			{
				if(rIM[t][i]<0 && IsResource(i)==false) //initial_marking[i]=0	 means P_A+end places. In the reversely firing, end places can be excluded.  D[i] are allowed to be 0. 					
				{
					p_pp=i;
					break; //Only one such place p_{i+1}.
				}
			}	
			if(writeFile==true)		
				outfile<<"p''=p_"<<p_pp+1<<endl;
			
				
			
			
			
			int maxRp_p, maxRp_pp;
			if(s_p.r1[p_p]>=s_p.r2[p_p])
				maxRp_p=s_p.r1[p_p];
			else
				maxRp_p=s_p.r2[p_p];
				
			if(s_p.r1[p_pp]>=s_p.r2[p_pp])
				maxRp_pp=s_p.r1[p_pp];
			else
				maxRp_pp=s_p.r2[p_pp];
				
				
				
			d=D[p_p]-maxRp_p;  //d is the time required for p_p to reversely enable t_{j+1} at s'. To reach the situation after 33333333333.
			if(writeFile==true)
				outfile<<"d="<<d<<endl;
				
			if(d>D[p_pp]-maxRp_pp)
			{	
				if(writeFile==true)
				{
					outfile<<"d>D[p_pp]-maxRp_pp!!!!!!!!!!!!!!!!!!!!!!"<<endl;
					outfile<<"--------------"<<endl;
				}
				continue;
			}
			
			
			if(D[p_p]==0 || maxRp_p>0 || p_p==p_pp)
				C &= Constraints(s_p, p_p, p_pp, d);
			else
			{
				bdd C1 = bddfalse;
				for(int k=0; k<=D[p_pp]-maxRp_pp-d; k++)
					C1 |= Constraints(s_p, p_p, p_pp, d+k);
				C &= C1;
			}
				
				
			//Try to find a state in Closed which satisfies C.
			bdd found = closed & C; 

			
			if(found != bddfalse)
			{
				//for debug
				if(writeFile==true)
				{
					outfile<<"p': p_"<<p_p+1<<endl;					
					outfile<<"d: "<<d<<endl;
					
					outfile<<"-------------------"<<endl;
					outfile<<"The number of satisfied state for C in CLOSED. bdd_satcount(found): "<<bdd_satcount(found)<<endl;
					outfile<<"They are: "<<fddset<<found<<endl;
					outfile<<"-------------------"<<endl;
				}
				
				bdd S_p=bdd_fullsatone(found);	//Note that S_p is a BDD. Not to use bdd_satone() since it may get more than one marking.	
				outfile<<"The number of state in S_p: "<<bdd_satcount(S_p)<<endl;
						
				bdd_allsat(S_p, BddToAStateHandler); //Transform S_p to sd.s by using BddToAStateHandler.	Uses abs()
				
				if(writeFile==true)
				{
					outfile<<"A satisfied state found in CLOSED:"<<endl;
					outfile<<"bdd_allsat(S_p, AllsatPrintHandler): ";					
					bdd_allsat(S_p, AllsatPrintHandler); //for debug
					outfile<<endl;
					
					outfile<<"outfile<<fddset<<S_p<<endl:";
					outfile<<fddset<<S_p<<endl<<endl;
				}
				
				
				int maxR, minR;				
				if(sd.s.mu[p_pp]==1)
				{
					if(sd.s.r1[p_pp]>=sd.s.r2[p_pp])	
						maxR=sd.s.r1[p_pp];
					else
						maxR=sd.s.r2[p_pp];
					sd.d=maxR;
				}
				else //sd.s.mu[p_pp]==2
				{
					if(sd.s.r1[p_pp]<=sd.s.r2[p_pp])	
						minR=sd.s.r1[p_pp];
					else
						minR=sd.s.r2[p_pp];
					sd.d=minR;
				}
				
				if(writeFile==true)
					outfile<<"The actural d is: "<<sd.d<<endl;
					
				return sd;				
			}//if(found != bddfalse)				
			
		}//if(rEt[j]==true)		
		
	}//for(int j=0;j<NT;j++)
	
	sd.s.f=-1; //f=-1 means no satisfied state in Closed.
	sd.d=0;  
	return sd;
	
}
			
			
			
//Constraints of F and R for the state s' such that s'[t>s whose time is d, \exists s'', s''[t>s', p' is the activity post-place of t', and p'' is the non-resource pre-place of t.
bdd Constraints(schedulingState s_p, int p_p, int p_pp, int d)			
{
			
			//Constraints for F.			
			int l=-1; //l is used to denote whether there is a loyal activity place which has tokens. If there is one (only one), l gives its index which starts from 0, otherwise, l=-1.
			
			
			//Cannot use for(int i=0;i<NP && s_p.mu[i]>0 && D[i]>0;i++) since i++ will only be executed after the body is performed.
			
				
			int N_loyal= sizeof(h3pj)/sizeof(h3pj[0]);  //obtain the number of elements in h3pj.
			for(int i=0;i<NP;i++)
			{
				if(s_p.mu[i]>0 && D[i]>0)
				{
					for(int k=0;k<N_loyal;k++)
						if(i+1==h3pj[k])  //The index of h3pj starts from 1.
							l=i; //at most one such a place
				}
			}
			
			if(writeFile==true)
			{
				outfile<<"l=p_"<<l+1<<endl;
				outfile<<"----------------"<<endl;
			}

			
			bdd C = bddtrue;
			if(l==-1)  //l==-1 means there is no such a place.
				C &= bvec_equ(F, bvec_con(MAXBITS,s_p.f-d));				
			else if(l==p_pp)  //if it is the output activity place of t'.
				C &= bvec_equ(F, bvec_con(MAXBITS,s_p.f));
			else //if it is an activity place other than the output one.
			{
				if(s_p.r1[l]>0 || s_p.r2[l]>0)
					C &= bvec_equ(F, bvec_con(MAXBITS,s_p.f));
				else //r1[l]=r2[l]=0
				{					
						int a=0;
						bdd C1=bddfalse;
						while(a<=d)
						{
							if(a==0)
								C1 |=  bvec_equ(F, bvec_con(MAXBITS,s_p.f-d)) & bvec_equ(R1[l], bvec_con(MAXBITS,0)) & bvec_equ(R2[l], bvec_con(MAXBITS,0));  //Note that for l
							else							
								C1 |=  bvec_equ(F, bvec_con(MAXBITS,s_p.f-d+a)) & ( bvec_equ(R1[l], bvec_con(MAXBITS,a)) | bvec_equ(R2[l], bvec_con(MAXBITS,a)));  
							a++;
						}
						C &= C1;					
					
				}
				
			}
				
			//Constraints of R in p', i.e., p_p.
			C &= bvec_equ(R1[p_p], bvec_con(MAXBITS,D[p_p])) | bvec_equ(R2[p_p], bvec_con(MAXBITS,D[p_p]));
			
			
			//Constraints of R in p'', i.e., p_pp.
			if(D[p_pp]>0)  //Since p'' may be a start place or an activity place with D=0.
			{
				if(s_p.mu[p_pp]==1)
				{
					
						if(d>0)
							C &= bvec_equ(R1[p_pp], bvec_con(MAXBITS,d)) | bvec_equ(R2[p_pp], bvec_con(MAXBITS,d));
						if(d==0)
							C &= bvec_equ(R1[p_pp], bvec_con(MAXBITS,0)) & bvec_equ(R2[p_pp], bvec_con(MAXBITS,0));
										
				}
				else  //s_p.mu[p_pp]==2
				{					
						C &= bvec_equ(R1[p_pp], bvec_con(MAXBITS,s_p.r1[p_pp]+d));
						C &= bvec_equ(R2[p_pp], bvec_con(MAXBITS,s_p.r2[p_pp]+d));
				}
			}
			
			
			//Constraints of R in other acitivity places.
			for(int i=0;i<NP;i++) 
			{
				if(s_p.mu[i]>0 && D[i]>0 && i!=p_p && i!=p_pp)
				{
					if(s_p.mu[i]==1)
					{
							if(s_p.r1[i]>0)
								C &= bvec_equ(R1[i], bvec_con(MAXBITS,s_p.r1[i]+d)) & bvec_equ(R2[i], bvec_con(MAXBITS,0));
							else if(s_p.r2[i]>0)
								C &= bvec_equ(R2[i], bvec_con(MAXBITS,s_p.r2[i]+d)) & bvec_equ(R1[i], bvec_con(MAXBITS,0));
							else
								C = C & ( bvec_lte(R1[i], bvec_con(MAXBITS,d)) & bvec_equ(R2[i], bvec_con(MAXBITS,0)) | bvec_lte(R2[i], bvec_con(MAXBITS,d)) & bvec_equ(R1[i], bvec_con(MAXBITS,0)));
					}
					else //s_p.mu[i]==2
					{
							if(s_p.r1[i]>0)
								C &= bvec_equ(R1[i], bvec_con(MAXBITS,s_p.r1[i]+d));
							else
								C &= bvec_lte(R1[i], bvec_con(MAXBITS,d));
							
							if(s_p.r2[i]>0)
								C &= bvec_equ(R2[i], bvec_con(MAXBITS,s_p.r2[i]+d));
							else
								C &= bvec_lte(R2[i], bvec_con(MAXBITS,d));
						
					}				
				}	
				
			}//for(int i=0;i<NP;i++) 
			
			
			return C;
}			



/********************************************************
// If you want to print the info into a outfile:

		ofstream outfile;
                outfile.open("result.txt");		//if want to append the file, then outfile.open("result.txt");
		Search(choice,0,0); //In Search(), use   outfile<<"Reachable Markings: "<<Petri<<endl;(cout<< is to display in screen)
		outfile.close();
**********************************************************/


int main()
{	
	cout<<"--------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"--------------------------Bo Huang, Nanjing University of Science and Technology, China-----------------------------"<<endl;
	cout<<"---------------------------------------------huangbo@njust.edu.cn---------------------------------------------------"<<endl;
	cout<<"--------------------------------------------------------------------------------------------------------------------"<<endl;
	cout<<"-----This program is to search a scheduling path from the initial state to a given final marking of a timed Petri"<<endl;
	cout<<"-----net of FMSs by using BDDs. Before running it, create a file of incidence matrix (netname_matrix.txt) and a file of initial marking"<<endl; 	
	cout<<"-----, processing time and final marking (netname_init.txt) in the same folder as the executable outfile. "<<endl; 
	cout<<"-----The parameters of the input net should be given and chosen at the begining of the code. "<<endl; 
	cout<<"-----It can search via BFS and the A* algorithm with a selected heuristic functions which should be specified by Search() in main()."<<endl;
	cout<<"--------------------------------------------------------------------------------------------------------------------"<<endl;

        cout<<"========================================================================================"<<endl;
	cout<<"ver1.1 (20170606_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Add comments;"<<endl;
	cout<<"2. Make the program more readable by changing some variable names and outputs;"<<endl;
	cout<<"3. Debug it for speed;"<<endl;	
	cout<<"4. When applied to a new PN, change the constant values of NP, NT, MAXNUMBER;"<<endl;
 	cout<<"   and MAXBITS at the beginning, need not to change the code;"<<endl;	
 	cout<<"----------------------------------------------------------------------------------------"<<endl;
 	cout<<"ver2.0 (20170622_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Revise for A* search with f=0;"<<endl;
	cout<<"2. Each element in (F,M,R1,R2) has the same number of bits;"<<endl;
	cout<<"3. The current version is for the number of tokens in any activity place is at most 2;"<<endl;
	cout<<"4. Rename all variables and functions;"<<endl;	
 	cout<<"----------------------------------------------------------------------------------------"<<endl;
 	cout<<"ver2.1 (20170707_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Define the constant int number_var_bigM =  3*NP+1;"<<endl;
	cout<<"2. It is for non-decreasing f;"<<endl;
	cout<<"3. Remove tfrom and Closed;"<<endl;
	cout<<"4. Add switch for debugging: const bool display_debug;"<<endl;
	cout<<"========================================================================================"<<endl;
	cout<<"ver2.2 (201707_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Add h, hmethod=0(h=0); hmethod=1(t_step-1, inadmissible);"<<endl;
	cout<<"2. Add BFS(breadth first search)!"<<endl;
	cout<<"3. Add hmethod=3(e_{Maxi}, admissible, consistent), it needs the given vector h3pj and h0!"<<endl;
	cout<<"4. Astar is for consistent h since f is considered to be non-decreasing!"<<endl;
	cout<<"========================================================================================"<<endl;
	cout<<"ver2.3 (201710_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Add CLOSED to compute the number of explored nodes and exclude the explored nodes from OPEN."<<endl;
	cout<<"2. Use bdd_printset() to replace bdd_allsat() in the display of results."<<endl;
	cout<<"========================================================================================"<<endl;
	cout<<"ver2.4 (20171121_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Add Construct() to find an optimal path in CLOSED and print it;"<<endl;
	cout<<"2. pro_time[] is replaced by D[];"<<endl;
	cout<<"========================================================================================"<<endl;
	cout<<"ver2.5 (20171207_Bo_Huang):"<<endl;
	cout<<"----------------------------------------------------------------------------------------"<<endl;
	cout<<"1. Revise Construct() to avoid deadend;"<<endl;
	cout<<"2. Rename variables according to the manuscript;"<<endl;
	cout<<"3. Output debug info to a outfile;"<<endl;
	cout<<"========================================================================================"<<endl;
	
	cout<<endl<<endl;
	
	struct timeval start,finish;  
	double totalTime;
	int optimal_f;
	
	writeFile=true;  //!!!!!!!!!!!!!!!
	
	
	if(writeFile == true)
	{
	        outfile.open("result.txt");	
	        outfile<<"Start searching..."<<endl<<endl;
	}
		
	

		cout<<"Start searching..."<<endl;
		InitMarking();
		
		//This program uses transposed incidence matrix.
		if(isIncidenceMatrix==true)
		{
			InitMatrix();	
			TransposeMatrix();
			//initCoInverseMatrix();
		}
		else
			InitTransposedMatrix();  
		
		
			
		
		//string choice = DisplayAndChoose();  
		
		//If you want to display the least info, use debug=false and least=true. 
		display_debug = false; //output the debug info.
		display_least = true; //output the least info.
		
		writeFile = true; //output the backtracking info to the outfile result.txt
		
		
		
		
		
	     	gettimeofday(&start,NULL);  				
	        
		//************************************************************************
		//Parameters of Search():
		//Parameter 1#(string, search method): "BFS"(breadth first search, NOT best first), "Astar"; 
		//Parameter 2#(int, hmethod): if 1#="BFS", then 0; if 1#="Astar", then: 0: h=0; 1: h=t_depth-1(only for post-t of operational places); 3: h=e_{Maxi};
		//Parameter 3#(int, h0): if 1#="BFS", then 0; if 1#="Astar": if hmethod=0, then h0=0; if hmethod =1, then h0=12 or else t_depth; if hmethod =2, then h0=0; if hmethod=3, then h0=h3_0, it also needs a predifined vector h3pj.
		//************************************************************************
		//Search("BFS",0, 0);//Use it and disable the rest if use BFS.
		Search("Astar",0, 0); //Use it and disable the rest if use h=0 in Astar.
		//optimal_f=Search("Astar",3, h3_0);  //Use it and disable the rest if use h=h_BDD in Astar. Note that only Astar return the optimal f value.	        
	        
	        gettimeofday(&finish,NULL);  
	     	totalTime = finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec) / 1000000.0;  
	     	
	     	
	  
	     	//cout<<"Your choice is: "<<choice<<endl;	
	     	
	     	//fdd_printset(open);
	     	//bdd_allsat(open, AllsatPrintHandler);
	     	cout<<endl;
	     	cout<<"-----------------------------------------"<<endl;
	     	cout<<"Search time: "<<totalTime<<" seconds."<<endl;       	
	     	cout<<"Scheduling states (F,M,R1,R2) # in CLOSED:"<<bdd_satcount(closed)<<endl;
	     	cout<<"Scheduling states (F,M,R1,R2) # in CLOSED (logarithm in base 2):"<<bdd_satcountln(closed)<<endl;
	     	cout<<"Node # of BDD for CLOSED:"<<bdd_nodecount(closed)<<endl<<endl; 	
	     	cout<<"Start backtracking..."<<endl;
	     	
	     	if(writeFile == true)
	     	{
		     	outfile<<endl;
		     	outfile<<"-----------------------------------------"<<endl;
		     	outfile<<"Search time: "<<totalTime<<" seconds."<<endl;       	
		     	outfile<<"Scheduling states (F,M,R1,R2) # in CLOSED:"<<bdd_satcount(closed)<<endl;
		     	outfile<<"Scheduling states (F,M,R1,R2) # in CLOSED (logarithm in base 2):"<<bdd_satcountln(closed)<<endl;
		     	outfile<<"Node # of BDD for CLOSED:"<<bdd_nodecount(closed)<<endl<<endl; 
		     	outfile<<"Start backtracking..."<<endl;	
	     	}
     	

     	
     	// Find the path in CLOSED
     	schedulingState s_0, s_G;
     	s_0.f=h3_0; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     	s_G.f=optimal_f; 

     	
     	for(int i=0;i<NP;i++)
     	{
     		s_0.mu[i]=initial_marking[i];
     		s_G.mu[i]=final_marking[i];
     		s_0.r1[i]=s_0.r2[i]=s_G.r1[i]=s_G.r2[i]=0;
     	}     	
     	
     	
     	gettimeofday(&start,NULL);        
        Construct(closed, s_0, s_G);
        gettimeofday(&finish,NULL);  
        
     	totalTime = finish.tv_sec - start.tv_sec + (finish.tv_usec - start.tv_usec) / 1000000.0;          
        cout<<"Backtracking time: "<<totalTime<<" seconds."<<endl;  
        
        if(writeFile == true)
        {		
		outfile<<"Backtracking time: "<<totalTime<<" seconds."<<endl;  
		outfile.close();
		
		cout<<"Result.txt has been generated!!!"<<endl;
	}
        
	return 0;	
}
