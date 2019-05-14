/**
 ** Definitions for structures 
 ** Statements of functions 
 ** @author:3160104633@zju.edu.com
 ** @date:2019/5/13
**/
#ifndef _DETECT_H_
#define _DETECT_H_

#define N 100
#define MAXM 5
#define MAXTIME 200
#define ABLE 1
#define DISABLED 0
#define LINE_MODEL 1
#define CIRCLE_MODEL 0
#define LENGTH(a) ((sizeof(a))/(sizeof(a[0])))

typedef struct _object Object;
typedef struct _object* ObjectPtr;
struct _object {
	int startcor[3]; //start cooridinate (x(t),y(t),z(t))
	int curcor[3]; //current cooridinate
	int oriseq; //object number, constant
	int seq; //object sequence, variable
	int velocity; 
	//for time denotation
	int begint;
	int endt;
	int sj;
	int pj;
	int second[2];
};

typedef struct _monitor Monitor;
typedef struct _monitor* MonitorPtr;
struct _monitor
{
	int center[3]; //center cooridinate (x(t),y(t),z(t))
	int id; //serial number for monitors
	int worklen;
	int slot[MAXTIME]; //work time for monitors
};

typedef struct _monitorset Mset;
struct _monitorset
{
	int amount;
	Monitor m[MAXM];
};

typedef struct _selected Selected;
typedef struct _selected* SelectedPtr;
struct _selected {
	int obseq;
	SelectedPtr next;
	SelectedPtr subs;//stores the sequence serial number for each item
	int subnum;
	int mark;
};

typedef struct _greedysolution Greedy;
typedef struct _greedysolution* GreedyPtr;
struct _greedysolution
{
	SelectedPtr objectarr;
	int capacity;
};

typedef struct _secondtracesolution SecondT;
struct _secondtracesolution
{
	int amount;
	int** secondset;
};

//void SecondTraceGreedy(int s[], int f[], int a[], int k);
SecondT* SecondTraceGreedy(ObjectPtr ob);

//void ExtendGreedySolution(GreedyPtr fs, ObjectPtr obarr,SecondT* sec);

//返回类型暂定为数组地址
int** getCoincidenceNum(SecondT* secres);

int** setCompatibleFactor(int** C, SecondT* secres);

double* setFactors(double percent, double thres, double* oriarr);

int* arrangeTarget(SecondT* secres, MonitorPtr m, int** comp, int am);

void readExcel(char* filename);

void writeExcel(int* res, int monitor);

void printResult(int* res,int monitor);

#endif // !_DETECT_H_
