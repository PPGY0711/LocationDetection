/**
** Implementation for functions
** @author:3160104633@zju.edu.com
** @date:2019/5/13
**/
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "detect.h"
#include "errormsg.h"
#define marked -1
#define NUM 1000 //设置0-1内随机数计算精度为小数点后4位
#define LENGTH(a) ((sizeof(a))/(sizeof(a[0])))

int aj[N], bj[N], pj[N], ts[N], te[N], v[N];
int seq_orinum[N][2];
int startP[N][3];

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
	int k,i;
	k = LENGTH(ob) - 1;
	for (i = 0; i < N; i++)
	{
		s[i] = ob[i].second[0];
		f[i] = ob[i].second[1];
	}
	greedy(s, f, res, &k);
	//the selected objects are marked by 1
	int* index;
	index = (int*)malloc(sizeof(int)*k);
	int j = 0;
	for (i = 0; i < N; i++) 
	{
		if (res[i] == 1)
		{
			index[j++] = i;
		}
	}

	//generate the first solution for greedy alogorithm
	SelectedPtr resarr;
	resarr = (SelectedPtr)malloc(sizeof(Selected)*k);
	memset(resarr, 0, sizeof(Selected)*k);
	
	for (i = 0; i < k-1; i++)
	{
		resarr[i].obseq = index[i];
		resarr[i].next = resarr+i+1;
	}
	resarr[k - 1].obseq = index[k - 1];
	resarr[k - 1].next = NULL;

	GreedyPtr fs = (GreedyPtr)malloc(sizeof(Greedy));
	fs->objectarr = resarr;
	fs->capacity = k;
	SecondT* sec;
	sec = ExtendGreedySolution(fs, s, f, index, k);
}

static void greedy(int s[], int f[], int a[], int* k)
{
	int i;
	int j = 0;
	a[0] = 1;
	int cnt = 1;
	for (i = 1; i < *k; i++)
	{
		if (s[i] > f[j])
		{
			a[i] = 1;
			j = i;
			cnt++;
		}
	}
	*k = cnt;
}

static SecondT* ExtendGreedySolution(GreedyPtr fs, int s[], int f[],int index[],int k)
{
	int i,j,id1,id2;
	int total = 1;
	for (i = 0, j = 1; j <= k - 1; i++, j++)
	{
		id1 = index[i];
		id2 = index[j];
		int subs;
		for (subs = id1; subs < id2; subs++)
		{
			//for the beginning of the loop, every item of selected sequence has the substitute of itself
			SelectedPtr tmp = fs->objectarr[i].subs;
			int* cnt = &(fs->objectarr[i].subnum);
			*cnt = 0;
			if (f[subs] < s[id2])
			{
				SelectedPtr node = new Selected;
				node->obseq = subs;
				node->next = node->subs = NULL;
				if (tmp == NULL) {
					tmp = node;
				}
				else {
					tmp->next = node;
				}
				tmp = tmp->next;
				*cnt = *cnt + 1;
			}
		}
		total = total*(fs->objectarr[i].subnum);
		if (j == k - 1)
		{
			for (subs = id2; subs < N; subs++)
			{
				SelectedPtr tmp = fs->objectarr[j].subs;
				int* cnt = &(fs->objectarr[j].subnum);
				*cnt = 0;
				if (f[id1] < s[subs])
				{
					SelectedPtr node = new Selected;
					node->obseq = subs;
					node->next = node->subs = NULL;
					if (tmp == NULL) {
						tmp = node;
					}
					else {
						tmp->next = node;
					}
					tmp = tmp->next;
					*cnt = *cnt + 1;
				}
			}
			total = total*(fs->objectarr[j].subnum);
		}
	}
	SecondT* resset = new SecondT;
	resset->amount = total;
	int** matrix = new int*[total];
	for (int n = 0; n < total; n++)
	{
		matrix[n] = new int[k];
		for (int m = 0; m < k; m++)
			matrix[n][m] = index[m];
	}
	
	//fs->objectarr[0].mark = marked;
	SelectedPtr tmp = fs->objectarr;
	for (int n = 0; n < total; n++)
	{
		int col = 0;
		int markcnt = 0;
		while (tmp&&markcnt<=n)
		{
			SelectedPtr cur = tmp->subs;
			//SelectedPtr nxt = tmp->next;
			while (cur)
			{
				if (cur->mark != marked)
				{
					matrix[n][col++] = cur->obseq;
					cur->mark = marked;
					markcnt++;
					break;
				}
				else
					markcnt++;
				cur = cur->next;
			}
			tmp = tmp->next;
		}
		tmp = fs->objectarr;
	}
	resset->secondset = matrix;
	return resset;
}

static int getIdealRes(int sortedaj[], int sortedbj[], int sortedpj[],int maxbj)
{
	//assume sum(bj[i]-aj[i]+pj[i] is sorted in ascending array
	int cnt = 0;
	int timesum = 0;
	while (timesum <= maxbj)
	{
		timesum += sortedbj[cnt] - sortedaj[cnt] + sortedpj[cnt];
		cnt++;
	}
	return --cnt;
}

//返回类型暂定为数组地址
int** getCoincidenceNum(SecondT* secres)
{
	int row, col;
	int** matrix = secres->secondset;
	row = LENGTH(matrix);
	col = LENGTH(matrix[0]);
	int** Cnum = new int*[row];
	for (int i = 0; i < row; i++)
	{
		Cnum[i] = new int[col];
		memset(Cnum[i], 0, sizeof(int)*col);
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			for (int k = j; k < col; k++)
			{
				//[ts,aj),[aj,bj]=>[ts,bj]
				if (((ts[matrix[i][j]] < ts[matrix[i][k]]) && (ts[matrix[i][k]] < bj[matrix[i][j]])) \
					|| ((ts[matrix[i][j]] < bj[matrix[i][k]]) && (bj[matrix[i][k]] < bj[matrix[i][j]])))
					Cnum[i][j] += 1;
			}
		}
	}
	return Cnum;
}

//对object对象进行基本数据处理
void ObjectHandle(ObjectPtr ob)
{
	handleObject(ob,aj,bj,pj,ts,te,v,startP,seq_orinum);
}

//目标数据处理内部接口
static void handleObject(ObjectPtr ob, int aj[], int bj[], int pj[],int ts[], int te[],int v[], int startP[][3], int seq_orinum[][2])
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

//计算可兼容性因子
int** setCompatibleFactor(int** C, SecondT* secres)
{
	int row, col;
	int** matrix = secres->secondset;
	row = LENGTH(matrix);
	col = LENGTH(matrix[0]);
	int** Comp = new int*[row];
	for (int i = 0; i < row; i++)
	{
		Comp[i] = new int[col];
		memset(Comp[i], 0, sizeof(int)*col);
	}
	srand(time(NULL));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			double weight[3] = { 0,0,0 };
			while (weight[2] <= 0) {
				for (int i = 0; i < 2; i++)
				{
					weight[i] = fabs(rand() % (NUM + 1) / (double)(NUM + 1));
				}
				weight[2] = 1 - weight[0] - weight[1];
			}
			//这个公式需要进一步完善
			Comp[i][j] = (weight[0] / C[i][j])*(weight[1] / pj[matrix[i][j]])*(weight[2] * (aj[matrix[i][j]] - pj[matrix[i][j]] - ts[matrix[i][j]]));
		}
	}
	return Comp;
}

//根据结果反馈调整权重因子
double* setFactors(double percent, double thres, double* oriarr)
{

}

int* arrangeTarget(SecondT* secres, MonitorPtr m, int** comp)
{
	//1.对任一解集，按照Comp值排序，Comp值大的目标优先放(排序后也要维护一个seq序列，否则就丢失了对应关系）
	//2.对于选中目标，将Sj设为能够放置的Sj最小值
	//3.第三块逻辑是查看当前需要放置的目标其区间是否已经被占用，若被占用，舍去，若未被占用，安排探测（已经写好）
	//23两点调用isAvailale()
}

//设置threshold
static double setThreshold(int ideal, int actual[]);

//资源是否被占用
static int isAvailable(MonitorPtr m, int obseq,int* _sj)
{
	int second = ABLE;
	int first = ABLE;
	for (int i = aj[obseq]; i <= bj[obseq]; i++)
	{
		if ((m->slot)[i] == DISABLED) {
			//second = DISABLED;
			return DISABLED;
		}
	}
	int i;
	for (i = ts[obseq]; i <= (aj[obseq] - pj[obseq]); i++)
	{
		int cnt = 0;
		int j = i;
		for ( ; j <= i + pj[obseq]; j++)
		{
			if(m->slot[j] == ABLE)
			cnt++;
		}
		if (cnt == pj[obseq]) {
			*_sj = j;
				return ABLE;
		}
	}
	return DISABLED;
}

//计算探测中心位置
void calculateCenterPos(int** targetList);

//读Excel文件接口
void readExcel(char* filename);

//写Excel工作表接口
void writeExcel(int* res, int monitor);

//将结果打印到终端
void printResult(int* res, int monitor);

//改写comp函数使sort从大到小排序
bool comp(const double &a, const double &b)
{
	return a>b;
}
