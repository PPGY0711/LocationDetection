/**
 ** Definitions for structures 
 ** Statements of functions 
 ** @name:detect.h
 ** @author:3160104633@zju.edu.com
 ** @date:2019/5/13
 **/
#ifndef _DETECT_H_
#define _DETECT_H_

#define N 100
#define MAXM 5
#define MAXTIME 200
#define ABLE 0
#define DISABLED -1
#define DEBUG 0
#define MANUAL 10
#define LINE_MODEL 1
#define CIRCLE_MODEL 0
#define LENGTH(a) ((sizeof(a))/(sizeof(a[0])))

typedef struct _object Object;
typedef struct _object* ObjectPtr;
struct _object {
	int startcor[3]; //start cooridinate (x(t),y(t),z(t))
//	int curcor[3]; //current cooridinate
	int oriseq; //object number, constant
	int seq; //object sequence, variable
	int velocity; 
	//for time denotation
	int begint;
	int endt;
	int sj;
	int pj;
	int second[2];
	int mark;
};

typedef struct _monitor Monitor;
typedef struct _monitor* MonitorPtr;
struct _monitor
{
	int center[3]; //center cooridinate (x(t),y(t),z(t))
	int id; //serial number for monitors
	int worklen;
	int position[MAXTIME][3];
	int* slot; //work time for monitors
};

typedef struct _monitorset Mset;
struct _monitorset
{
	int amount;
	Monitor m[MAXM];
};

//被选中的目标，subs存放其分区目标 
typedef struct _selected Selected;
typedef struct _selected* SelectedPtr;
struct _selected {
	int obseq;
	SelectedPtr next;
	SelectedPtr subs;
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

//二次追踪贪心算法解集 
typedef struct _secondtracesolution SecondT;
struct _secondtracesolution
{
	int amount;
	int* secondset;
};


//贪心算法内部接口 
static int greedy(int s[], int f[], int a[], int k);
//扩展解集函数（内部接口） 
static SecondT* ExtendGreedySolution(GreedyPtr fs, int s[], int f[], int index[]);
//获取理想结果数 
static void getIdealRes(int totaltime[], int maxbj, int* MAXN); 
static void MaxArrangeNumber(int aj[], int bj[], int pj[], int* MAXN);
//目标数据处理内部接口
static void handleObject(ObjectPtr ob, int aj[], int bj[], int pj[], int ts[], int te[], int v[], int startP[][3], int seq_orinum[][2]);
//探测器内容初始化 
static void initialMonitors(MonitorPtr m,int id);
//获取当前未被追踪的目标集合 
static SecondT* getUntracedSet();
//资源是否被占用
static int isAvailable(MonitorPtr m, int obseq, int* _sj,int flag);
//获取当前状态下所有探测器的剩余可用时间区间 
static void getAvailableSlot(Mset* ms);
//检查剩余区间Rem是否可以成功追踪同伴 
static void CheckRem(MonitorPtr m,SecondT* secres);
//打乱可替换项顺序 
static void shuffle(GreedyPtr aim);
//计算成功追踪目标数 
static int getTracedNum(int* arr); 
//打乱顺序函数 
static int* random(int n,int i);
//计算探测中心位置
void calculateCenterPos(int* targetList,MonitorPtr m,int AN,int id);
//按照二次追踪时间实现的贪心算法 
SecondT* SecondTraceGreedy(ObjectPtr ob);
//对object对象进行基本数据处理
void ObjectHandle(ObjectPtr ob);
//得到所有探测器的安排
void getArrangement(Mset* monitors, SecondT* secres, int amount, int flag);
//对于同一个探测器安排其调度 
int* arrangeTarget(SecondT* secres, MonitorPtr m,int* AN,int flag);
//读Excel文件接口
void readExcel();
//打印结果 
void printResult(int* resarr, int id,int AN,int flag);
bool compare(const int &a, const int &b);
bool compare1(const int &a, const int &b);

#endif // !_DETECT_H_
