#pragma warning(disable:4996)
/**
 ** Implementation for functions
 ** @name:detect.cpp
 ** @author:3160104633@zju.edu.com
 ** @date:2019/5/13
 **/
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <string.h>
#include <streambuf>
#include <fstream>
#include <ctime>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "detect.h"
#include "errormsg.h"
using namespace std;

#define PI 3.14
#define TRYTIME 1
#define marked -1
#define NOOBJECT -1
#define BEENTRACED -1
#define NUM 1000 //设置0-1内随机数计算精度为小数点后4位
#define LENGTH(a) ((sizeof(a))/(sizeof(a[0])))

//通过目标数据可以得到的基本数据
Object obarr[100],tmpobarr[100];
double weight[3] = { 0,0,0 };
int aj[N], bj[N], pj[N], ts[N], te[N], v[N];
int seq_orinum[N][2];
int startP[N][3];
GreedyPtr fs,remfs;
int colnum;
int* mid;
int sjSet[6 ][N];
int ASlot[MAXTIME];

//按照二次追踪时间实现的贪心算法 
SecondT* SecondTraceGreedy(ObjectPtr ob)
{
	//the object array has already been sorted according to end time of second trace
	int* s;
	int* f;
	int* res;
	s = (int*)malloc(sizeof(int)*N);
	f = (int*)malloc(sizeof(int)*N);
	res = (int*)malloc(sizeof(int)*N);
	if (s == NULL || f == NULL || res == NULL) {
		errorReport("Can't mallocating room for array.");
		exit(-1);
	}
	memset(res, 0, sizeof(int)*N);
	int k, i;
	k = LENGTH(obarr);
	for (i = 0; i < N; i++)
	{
		s[i] = ob[i].second[0];
		f[i] = ob[i].second[1];
	}
	k = greedy(s, f, res, k);
	colnum = k;
	int* index;
	index = (int*)malloc(sizeof(int)*k);
	int j = 0;
	for (i = 0; i < N; i++)
	{
		if (res[i] == 1)
		{
			index[j] = i;
			j++;
		}
	}

	//generate the first solution for greedy alogorithm
	SelectedPtr resarr;
	resarr = (SelectedPtr)malloc(sizeof(Selected)*k);
	memset(resarr, 0, sizeof(Selected)*k);

	for (i = 0; i < k - 1; i++)
	{
		resarr[i].obseq = index[i];
		resarr[i].next = &resarr[i + 1];
	}
	resarr[k - 1].obseq = index[k - 1];
	resarr[k - 1].next = NULL;
	SelectedPtr tmp = resarr;
	fs = (GreedyPtr)malloc(sizeof(Greedy));
	fs->objectarr = resarr;
	fs->capacity = k;
	SecondT* sec = new SecondT;
	//扩展二次追踪最优解集

	sec = ExtendGreedySolution(fs, s, f, index);
	return sec;
}

static int greedy(int s[], int f[], int a[], int k)
{
	int i,z;
	for(i = 0; i < N;i++)
	{
		if(tmpobarr[i].mark!=BEENTRACED){
			a[i] = 1;
			break;
		}
			
	}
	int cnt = 1;
	for (z = i+1; z < k; z++)
	{
		if (s[z] > f[i] && tmpobarr[z].mark!=BEENTRACED)
		{
			a[z] = 1;
			i = z;
			cnt++;
		}
	}
	return cnt;
}

//扩展解集函数 
static SecondT* ExtendGreedySolution(GreedyPtr fs, int s[], int f[], int index[])
{
	int i, j, id1, id2;
	
	for (i = 0, j = 1; j <= colnum - 1; i++, j++)
	{
		id1 = index[i];
		id2 = index[j];
		int subs;
		int* cnt = &(fs->objectarr[i].subnum);
		*cnt = 0;
		SelectedPtr tmp = fs->objectarr[i].subs;
		for (subs = id1; subs < id2; subs++)
		{
			//for the beginning of the loop, every item of selected sequence has the substitute of itself
			//solution B) term:
			if (f[subs] < s[id2] && s[subs] >= s[id1])
			//solution A) term:
			//if(f[subs] < s[id2] )
			{
				SelectedPtr node = new Selected;
				node->obseq = subs;
				node->next = node->subs = NULL;
				if (tmp == NULL) {
					fs->objectarr[i].subs = node;
					tmp = fs->objectarr[i].subs;
				}
				else {
					tmp->next = node;
					tmp = tmp->next;
				}
				node->mark = *cnt;
				*cnt = *cnt + 1;
			}
		}
		if (j == colnum - 1)
		{
			int* cnt = &(fs->objectarr[j].subnum);
			*cnt = 0;
			SelectedPtr tmp = fs->objectarr[j].subs;
			for (subs = id2; subs < N; subs++)
			{
				//solution B) term:
				if (f[id1] < s[subs] && s[subs] >= s[id1])
				//solution A) term:
				//if(f[id1] < s[subs])
				{
					SelectedPtr node = new Selected;
					node->obseq = subs;
					node->next = node->subs = NULL;
					if (tmp == NULL) {
						fs->objectarr[j].subs = node;
						tmp = fs->objectarr[j].subs;
					}
					else {
						tmp->next = node;
						tmp = tmp->next;
					}
					node->mark = *cnt;
					*cnt = *cnt + 1;
				}
			}
		}
	}

	SecondT* resset = new SecondT;
	resset->amount = TRYTIME;
	int* arr = new int[colnum];
	for (int n = 0; n < colnum; n++)
	{
		arr[n] = index[n];
	}
	Selected tmp = fs->objectarr[0];
	int col = 0;
	while ((tmp.next) != NULL)
	{
		SelectedPtr cur = tmp.subs;	
		int choose = 1;
		for (int cnt = 1; cnt < choose; cnt++)
		{
			cur = cur->next;
		}
		arr[col] = cur->obseq;
		col++;
		SelectedPtr t = tmp.next;
		tmp = *t;
	}
	resset->secondset = arr;
	return resset;
}

static void getIdealRes(int totaltime[], int maxbj, int* MAXN)
{
	//assume sum(bj[i]-aj[i]+pj[i] is sorted in ascending array
	int cnt = 0;
	int timesum = 0;
	while (timesum <= maxbj)
	{
		timesum += totaltime[cnt];
		cnt++;
	}
	*MAXN = --cnt;
}

static void MaxArrangeNumber(int aj[], int bj[], int pj[], int* MAXN)
{
	int* maxbj = max_element(bj, bj + N);
	int totaltime[N];
	for (int i = 0; i < N; i++)
	{
		if(tmpobarr[i].mark != BEENTRACED)
			totaltime[i] = bj[i] - aj[i] + pj[i];
	}
	sort(totaltime, totaltime + N, compare1);
	getIdealRes(totaltime, *maxbj, MAXN);
}

//对object对象进行基本数据处理
void ObjectHandle(ObjectPtr ob)
{
	handleObject(ob, aj, bj, pj, ts, te, v, startP, seq_orinum);
}

//目标数据处理内部接口
static void handleObject(ObjectPtr ob, int aj[], int bj[], int pj[], int ts[], int te[], int v[], int startP[][3], int seq_orinum[][2])
{
	for (int i = 0; i < N; i++)
	{
		aj[i] = ob[i].second[0];
		bj[i] = ob[i].second[1];
		pj[i] = ob[i].pj;
		ts[i] = ob[i].begint;
		te[i] = ob[i].endt;
		v[i] = ob[i].velocity;
		for (int j = 0; j < 3; j++)
		{
			startP[i][j] = ob[i].startcor[j];
		}
		seq_orinum[i][0] = ob[i].seq;
		seq_orinum[i][1] = ob[i].oriseq;
	}
}

//探测器内容初始化 
static void initialMonitors(MonitorPtr m,int id)
{
	m->id = id;
	m->slot = new int[MAXTIME];
	for(int i = 0; i < MAXTIME; i++)
	{
		m->slot[i] = 0;
	}
	for(int i = 0; i < MAXTIME; i++)
	{
		m->position[i][0] = m->position[i][1] =m->position[i][2] =0;
	}
}

//获取当前未被追踪的目标集合 
static SecondT* getUntracedSet()
{
	return SecondTraceGreedy(tmpobarr);
}

//得到所有探测器的安排
void getArrangement(Mset* monitors, SecondT* secres, int amount, int flag)
{
	
	for (int i = 0; i < amount; i++)
	{
		initialMonitors(&monitors->m[i],i+1);
		int AN = 0;
		int* resarr = NULL;
		if(i != 0) 
			secres = getUntracedSet();
		resarr = arrangeTarget(secres, &(monitors->m[i]), &AN,flag);
		//calculate coordinate	
		calculateCenterPos(resarr,&monitors->m[i],AN,i+1); 
		//print result
		printResult(resarr,i+1,AN,0);
	}
}

//对于同一个探测器安排其调度 
int* arrangeTarget(SecondT* secres, MonitorPtr m,int* AN,int flag)
{
	//对当前的探测器m进行操作,secres选出来的是我当前选择的一组二次不互相干扰的目标序列 
	int MAXN = 0;//考虑所有的目标，按照总时间长短得到的能探测的目标的理论最大值
	MaxArrangeNumber(aj, bj, pj, &MAXN);
	//初始化结果列表
	int* resarr;
	resarr = new int[MAXN];
	memset(resarr, NOOBJECT, sizeof(int)*MAXN);
	//初始化sj列表
	int curcol = 0;
	int* tmp = secres->secondset;
	for (int k = 0; k < colnum; k++)
	{
		//对一个探测器进行安排
		SelectedPtr cur = fs->objectarr[k].subs;
		SelectedPtr ppos = cur;
		SelectedPtr precur = NULL;
		SelectedPtr pos = NULL; 
		if(k > 0){
			precur = fs->objectarr[k-1].subs;
			pos = precur;
		}
		while(cur != NULL)
		{	
			//补偿措施：当当前的目标不能放进去时，在其替换项中按顺序找出可替换项安排进去
			if(tmpobarr[cur->obseq].mark != BEENTRACED && isAvailable(m,cur->obseq,&sjSet[m->id-1][cur->obseq],flag) == ABLE)
			{
				tmpobarr[cur->obseq].mark = BEENTRACED;
				resarr[curcol] = cur->obseq;
				curcol++;
				break;
			}
			cur = cur->next;
		}
		cur = ppos;
	}
	*AN = curcol;
	return resarr;
}

//资源是否被占用
static int isAvailable(MonitorPtr m, int obseq, int* _sj,int flag)
{
	if(flag == 0)
	{	
		for (int i = aj[obseq]; i <= bj[obseq]; i++)
		{
			if (m->slot[i] == DISABLED) {
				return DISABLED;
			}
		}
		for (int i = ts[obseq]; i <= (aj[obseq] - pj[obseq]); i++)
		{
			int cnt = 0;
			int j = i;
			//选sj的位置 
			for (; j < i + pj[obseq]; j++)
			{
				if (m->slot[j] == ABLE)
					cnt++;
			}
			if (cnt == pj[obseq]) {
				*_sj = i;
				for (int z = aj[obseq]; z <= bj[obseq]; z++)
					m->slot[z] = DISABLED;
				for (int z = *_sj; z < *_sj + pj[obseq]; z++)
					m->slot[z] = DISABLED;
				return ABLE;
			}
		}
		return DISABLED;
	}
	else
	{
		for (int i = aj[obseq]; i <= bj[obseq]; i++)
		{
			if (m->slot[i] == DISABLED) {
				return DISABLED;
			}
		}
		for (int i = ts[obseq]; i <= (aj[obseq] - pj[obseq]); i++)
		{
			int cnt = 0;
			mid = new int[pj[obseq]+bj[obseq]-aj[obseq]+1];
			int j = i;
			//选sj的位置 
			for (; j < i + pj[obseq]; j++)
			{
				if (m->slot[j] >= 0){
					//mid记录当前时刻由探测器m[id]探测该目标 
					mid[cnt] = m->slot[j]+1;
					cnt++;
				}
			}
			if (cnt == pj[obseq]) {
				*_sj = i;
				for (int z = aj[obseq]; z <= bj[obseq]; z++)
				{
					mid[pj[obseq]+z-aj[obseq]] = m->slot[z]+1;
					m->slot[z] -= 1;
				}
				for (int z = *_sj; z < *_sj + pj[obseq]; z++)
					m->slot[z] -= 1;
				return ABLE;
			}
		}
		return DISABLED;
	}	
}

//获取当前状态下所有探测器的剩余可用时间区间 
static void getAvailableSlot(Mset* ms)
{
	
	for(int i = 0; i < MAXTIME;i++)
	{
		ASlot[i] = -1;
	}
	for(int i = 0;i < ms->amount;i++)
	{
		for(int j = 0; j < MAXTIME;j++)
		{
			if(ms->m[i].slot[j]==0)
				ASlot[j]++;
		}
	}
}

//检查剩余区间Rem是否可以成功追踪同伴 
static void CheckRem(MonitorPtr m,SecondT* secres)
{
	int AN = 0;
	int* resarr = NULL;
	secres = getUntracedSet();
	resarr = arrangeTarget(secres, m, &AN,1);
	//calculate coordinate	
	calculateCenterPos(resarr,m,AN,6);
	//print result
	printResult(resarr,m->id,AN,1);
	
}


//打乱可替换项顺序 
static void shuffle(GreedyPtr aim)
{
	GreedyPtr tmp = aim;
	for(int i = 0;i < colnum;i++)
	{
		int* newind = random(tmp->objectarr[i].subnum+1,i);
		SelectedPtr head = new Selected;
		head->next = NULL;
		SelectedPtr ptr = head;
		SelectedPtr ppos = tmp->objectarr[i].subs;
		for(int j = 0; j < tmp->objectarr[i].subnum; j++ )
		{
			SelectedPtr pos = tmp->objectarr[i].subs;
			SelectedPtr node = new Selected;
			while(pos != NULL){
					
				if(pos->mark == newind[j] && head->next != NULL)
				{
					node->obseq = pos->obseq;
					node->mark = pos->mark;
					if(newind[j] == 0)
						node->subnum = pos->subnum;
					ptr->next = node;
					ptr = node;
				}
				if(pos->mark == newind[j] && head->next == NULL){
					node->obseq = pos->obseq;
					node->mark = pos->mark;
					if(newind[j] == 0)
						node->subnum = pos->subnum;
					head->next = node;
					ptr = node;
				}
				pos = pos->next;
			}
		}
		ppos = head->next;
	}
} 

//计算成功追踪目标数 
static int getTracedNum(int* arr)
{
	int cnt = 0;
	for (int i = 0; i < colnum; i++)
	{
		if (arr[i] == NOOBJECT)
			break;
		cnt++;
	}
	return cnt;
}

//计算探测中心位置
void calculateCenterPos(int* targetList,MonitorPtr m,int AN,int id)
{
	int time = 0; 
	for(int i = 0;i < AN; i++)
	{
		int tv = v[targetList[i]];
		int taj = aj[targetList[i]];
		int tbj = bj[targetList[i]];
		int tts = ts[targetList[i]];
		int tpj = pj[targetList[i]];
		int tsj = sjSet[id-1][targetList[i]];
#if LINE_MODEL		
		//目标坐标 
		int tmppos[3] ={0};
		//探测器坐标 
		//由于物体运动速度≤10/s，每秒改变一次探测器位置可以保证完成追踪 

		for(int j = 0; j < 3;j++)
		{
			if(j == 0){
				tmppos[j] = startP[targetList[i]][j]+(tsj-tts)*tv;
				m->position[tsj][j] = tmppos[j]+5;
			}
			tmppos[j] = startP[targetList[i]][j];
			m->position[tsj][j] = tmppos[j];
		}
		for(time = tsj; time < tsj+tpj; time++)
		{
			tmppos[0] = tmppos[0]+(time-tsj)*tv;
			m->position[time][0] = tmppos[0]+5;
		}
#endif
#if CIRCLE_MODEL
		//设所有的目标从同一直线的不同高度抛出，探测中心位置就是物体位置，不再以秒为单位变化而是时刻变化，这里以秒为单位记录
		//返回三维柱坐标（初始高度，开始时间，初速度，速度&z轴夹角，速度极坐标方向，当前时间） (r,phi,z)->(0,1,2)
		
		srand(time(NULL));
		double theta = (rand()%(PI))+1;
		double tmp1 = v[targetList[i]]*sin(theta)*(tsj-tts)-0.5*9.8*(tsj-tts)*(tsj-tts);
		double tmp2 = tmp1-9.8*(tsj-tts)-0.5*9.8;
		m->position[tsj][2] = (tmp1+tmp2)/2;
		tmp1 = v[targetList[i]]*cos(theta)*(tsj-tts);
		tmp2 = tmp1+v*cos(theta);
		m->position[tsj][0] = (tmp1+tmp2)/2;
		for(time = tsj; time < tsj+tpj; time++)
		{
			theta = (rand()%(PI))+1;
			tmp1 = v[targetList[i]]*sin(theta)*(time-tsj)-0.5*9.8*(time-tsj)*(time-tsj);
			tmp2 = tmp1-9.8*(time-tsj)-0.5*9.8;
			m->position[time][2] = (tmp1+tmp2)/2;
			tmp1 = v[targetList[i]]*cos(theta)*(time-tsj);
			tmp2 = tmp1+v*cos(theta);
			m->position[time][0] = (tmp1+tmp2)/2;
		}
#endif
	}
}

//读Excel文件接口
void readExcel()
{
#if LINE_MODEL
	char* data[8];
#endif
#if CIRCLE_MODEL
	char* data[10];
#endif
	ifstream iFile("Problem A.csv");
	if (!iFile.is_open())
	{
		errorReport("Can't open .csv file.");
		exit(-1);
	}
	char buffer[256];
	int i = 0, j = 0;
	char* tmp;
	const char *d = ",";
	while (!iFile.eof())
	{
		iFile.getline(buffer, 100);
		tmp = strtok(buffer, d);
		if (1 <= i && i <= 100) {
			while (tmp)
			{
				data[j++] = tmp;
				tmp = strtok(NULL, d);
			}
			obarr[i - 1].mark = 0;
			obarr[i - 1].seq = i - 1;
			obarr[i - 1].oriseq = atoi(data[0]);
			obarr[i - 1].begint = atoi(data[1]);
			obarr[i - 1].endt = atoi(data[2]);
			obarr[i - 1].startcor[0] = atoi(data[3]);
#if LINE_MODEL
			obarr[i - 1].startcor[1] = obarr[i - 1].startcor[2] = 0;
			obarr[i - 1].velocity = atoi(data[4]);
			obarr[i - 1].pj = atoi(data[5]);
			obarr[i - 1].second[0] = atoi(data[6]);
			obarr[i - 1].second[1] = atoi(data[7]);
#endif
#if CIRCLE_MODEL
			obarr[i - 1].startcor[1] = atoi(data[4]);
			obarr[i - 1].startcor[2] = atoi(data[5]);
			obarr[i - 1].velocity = atoi(data[6]);
			obarr[i - 1].pj = atoi(data[7]);
			obarr[i - 1].second[0] = atoi(data[8]);
			obarr[i - 1].second[1] = atoi(data[9]);
#endif	
		}
		j = 0;
		i++;
	}
}

//打印结果 
void printResult(int* resarr, int id,int AN,int flag)
{
	if(flag == 0)
	{
		printf("Monitor [%d] traced Actual num %d target:\n" ,id ,AN);
		printf("Target id\tTarget seq\tsj\tsj+pj\taj\tbj\n");
		for(int j = 0;j < AN; j++)
		{
			printf("\t%d\t\t%d\t%d\t%d\t%d\t%d\n",obarr[resarr[j]].oriseq,resarr[j],sjSet[id-1][resarr[j]],sjSet[id-1][resarr[j]]+pj[resarr[j]]-1,aj[resarr[j]],bj[resarr[j]]);
		}
	}
	else
	{
		
		printf("Monitor new traced Actual num %d target:\n" ,AN);
		printf("Target id\tTarget seq\tsj\tsj+pj\taj\tbj\n");
		for(int j = 0;j < AN; j++)
		{
			printf("\t%d\t\t%d\t%d\t%d\t%d\t%d\n",obarr[resarr[j]].oriseq,resarr[j],sjSet[id-1][resarr[j]],sjSet[id-1][resarr[j]]+pj[resarr[j]]-1,aj[resarr[j]],bj[resarr[j]]);
			
			for(int k = 0; k < pj[resarr[j]]+bj[resarr[j]]-aj[resarr[j]]+1;k++)
			{
				if(k==0)
					printf("In [sj,sj+pj]:\n");
				if(k==pj[resarr[j]])
					printf("In [aj,bj]:\n");
				printf("This moment traced by monitor[%d] \n",mid[k]);

			}
			
			printf("\n");
		}
	}
}


bool compare(const int &a, const int &b)
{
	return a>b;
}

bool compare1(const int &a, const int &b)
{
	return a<b;
}
 
static int* random(int n,int i)
{
	
	int index, tmp, k;
	srand(time(NULL));
	int* a = new int[n-1];
	int* b = new int[n-1];
	SelectedPtr node = fs->objectarr[i].subs;
	for(k = 0;k < n-1;k++)
	{
		a[k] = bj[node->obseq]-aj[node->obseq]+pj[node->obseq];
		b[k] = a[k];
		node = node->next;
	}
	sort(a,a+n-1,compare);
	int* c = new int[n-1];
	int j = 0;
	for(i = 0;i< n-1;i++)
	{
		for(j = 0;j < n-1;j++)
		{
			if(b[i] == a[j])
				c[i] = j;
		}
	}
	return c;
}

int main()
{
	//探测器数量
	int amount;
	//探测器集合
	Mset* m = new Mset;
	//二次追踪贪心算法解集
	SecondT* secset = NULL;

	readExcel();
	//times控制循环次数 
	
	int times;
	int i = 0;
	if(DEBUG) 
		times = 1;
	else
	{
		printf("==Please input loop times: ");
		scanf("%d",&times);
	}
	while (i<times) {
		for (int i = 0; i < N; i++)
		{
		tmpobarr[i] = obarr[i];
		}
		ObjectHandle(tmpobarr);
		if (DEBUG)
		{
			printf("==Monitor amount: 5\n");
			m->amount = 5;
			amount = 5;
		}
		else
		{
			printf("==Please input Monitor amount: ");
			scanf("%d", &amount);
			m->amount = amount;	
		}
		secset = SecondTraceGreedy(tmpobarr);
		shuffle(fs);
		getArrangement(m, secset, amount, 0);
		if(amount == 5){
			//收集剩余区间 
			Monitor tmpm;
			getAvailableSlot(m);
			tmpm.id = 6;
			tmpm.slot = ASlot;
			CheckRem(&tmpm,secset);		
		}
		i++;
	}
	
	system("PAUSE");
	return 0;
}
