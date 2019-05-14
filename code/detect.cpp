/**
** Implementation for functions
** @author:3160104633@zju.edu.com
** @date:2019/5/13
**/
#include <iostream>
#include <stdlib.h>
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
#define marked -1
#define NOOBJECT -1
#define BEENTRACED -1
#define NUM 1000 //设置0-1内随机数计算精度为小数点后4位
#define LENGTH(a) ((sizeof(a))/(sizeof(a[0])))

//通过目标数据可以得到的基本数据
Object obarr[100];
double weight[3] = { 0,0,0 };
int aj[N], bj[N], pj[N], ts[N], te[N], v[N];
int seq_orinum[N][2];
int startP[N][3];
int** sjSet;

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

static void getIdealRes(int totaltime[],int maxbj,int* MAXN)
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

static void MaxArrangeNumber(int aj[], int bj[], int pj[],int* MAXN)
{
	int* maxbj;
	*maxbj = *max_element(bj, bj + N);
	int totaltime[N];
	for(int i = 0; i < N; i++)
	{
		totaltime[N] = bj[i] - aj[i] + pj[i];
	}
	sort(totaltime, totaltime + N, compare);
	getIdealRes(totaltime, *maxbj, MAXN);
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
				if (matrix[i][j]!= BEENTRACED && \
					((ts[matrix[i][j]] < ts[matrix[i][k]]) && (ts[matrix[i][k]] < bj[matrix[i][j]])) \
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
		//sj[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			startP[i][j] = ob[i].startcor[j];
		}
		seq_orinum[i][0] = ob[i].seq;
		seq_orinum[i][1] = ob[i].oriseq;
	}
}

//计算可兼容性因子
int** setCompatibleFactor(int** C, SecondT* secres, int state, int CHANGE)
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
	if (state == 0 && CHANGE ==0) {
		while (weight[2] <= 0) {
			for (int i = 0; i < 2; i++)
			{
				weight[i] = fabs(rand() % (NUM + 1) / (double)(NUM + 1));
			}
			weight[2] = 1 - weight[0] - weight[1];
		}
	}
	if (CHANGE)
	{
		//set new weights on the basis of the original values
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			//这个公式需要进一步完善
			if(matrix[i][j] != BEENTRACED)
				Comp[i][j] = (weight[0] / C[i][j])*(weight[1] / pj[matrix[i][j]])*(weight[2] * (aj[matrix[i][j]] - pj[matrix[i][j]] - ts[matrix[i][j]]));
		}
	}
	return Comp;
}

//根据结果反馈调整权重因子
double* setFactors(double percent, double thres, double* oriarr)
{
	//还没实现
	return oriarr;
}

//得到所有探测器的安排
void getArrangement(Mset monitors, SecondT* secres,int amount,int CHANGE)
{
	//初始化探测器数据
	monitors.amount = amount;
	for (int i = 0; i < amount; i++)
	{
		monitors.m[i].id = i;
		monitors.m[i].slot = new int*[LENGTH(secres->secondset)];
		for (int j = 0; j < LENGTH(secres->secondset); j++)
		{
			monitors.m[i].slot[j] = new int[MAXTIME];
			memset(monitors.m[i].slot[j], ABLE, sizeof(int)*MAXTIME);
		}
		//重新计算Comp指标，重合数，但不改变随机权重的值，一组随机数算一整组，之后的权重是根据结果好坏来调整的
		int** C = getCoincidenceNum(secres);
		if (CHANGE)
			refreshMonitors(monitors);
		int** Comp = setCompatibleFactor(C, secres,i,CHANGE);
		int** sorted = sortedCompandSeq(secres->secondset, Comp);
		int** resmat = arrangeTarget(secres, &(monitors.m[i]), sorted);
		int* index;
		int scnt;
		int** finalRes = findBestSolution(resmat,index,&scnt);
		int** finalsjSet = getsjSet(sjSet,index,scnt);
		//int** TPosition = calculateCenterPos(finalRes,v,aj,bj,startP, finalsjSet);
		//printResult()
	}
}

static int** getsjSet(int** sjSet, int index[],int scnt)
{
	int no = scnt;
	int** finalsjSet = new int*[scnt];
	for (int i = 0; i < scnt; i++)
	{
		finalsjSet[i] = sjSet[index[i]];
	}
	return finalsjSet;
}

static void refreshMonitors(Mset monitors)
{
	for (int i = 0; i < monitors.amount; i++)
	{
		//monitors.m[i].id = i;
		//monitors.m[i].slot = new int*[LENGTH(secres->secondset)];
		for (int j = 0; j < MAXTIME; j++)
		{
			//monitors.m[i].slot[j] = new int[MAXTIME];
			memset(monitors.m[i].slot[j], ABLE, sizeof(int)*MAXTIME);
		}
	}
}

int** arrangeTarget(SecondT* secres, MonitorPtr m, int** CsortSeq)
{

	int MAXN;//考虑所有的目标，按照总时间长短得到的能探测的目标的理论最大值
	MaxArrangeNumber(aj, bj, pj, &MAXN);
	//初始化结果矩阵
	int** resmat = new int*[LENGTH(secres->secondset)];
	for (int j = 0; j < LENGTH(secres->secondset); j++)
	{
		resmat[j] = new int[MAXN];
		memset(resmat[j], NOOBJECT, sizeof(int)*MAXN);
	}
	//初始化sj矩阵
	int** sj = new int*[LENGTH(secres->secondset)];
	for (int j = 0; j < LENGTH(secres->secondset); j++)
	{
		sj[j] = new int[LENGTH(secres->secondset[0])];
		memset(resmat[j], 0, sizeof(int)*LENGTH(secres->secondset[0]));
	}
	//对于选中目标，将Sj设为能够放置的Sj最小值
	for (int j = 0; j < LENGTH(secres->secondset); j++)
	{
		int curcol = 0;
		for (int k = 0; k < LENGTH(secres->secondset[0]); k++)
		{
			if (CsortSeq[j][k] != BEENTRACED && isAvailable(m, CsortSeq[j][k], &sj[j][CsortSeq[j][k]],j) == ABLE)
			{
				resmat[j][curcol++] = CsortSeq[j][k];
				//1.从集合中去掉已经排到前一个monitor的目标
				//在secondset里面把seq等于CsortSeq的seq置为-1
				for (int z = 0; z < LENGTH(secres->secondset[0]); z++)
				{
					if (secres->secondset[j][z] == CsortSeq[j][k])
						secres->secondset[j][z] = BEENTRACED;
				}
				CsortSeq[j][k] = BEENTRACED;
			}
		}
		//对每个解集重置slot
		//memset(m->slot, ABLE, sizeof(int)*MAXTIME);
	}
	sjSet = sj;
	return resmat;
}

static int** sortedCompandSeq(int** secondset, int** Comp)
{
	//返回一个二维数组，行数为secondset行数，列数为secondset列数，但是记录的是按照Comp值从大到小sort之后的seq值。
	int** Csort = new int*[LENGTH(secondset)];
	for (int i = 0; i < LENGTH(secondset); i++)
	{
		Csort[i] = new int[LENGTH(secondset[0])];
		//memset(Csort[i], 0, sizeof(int)*LENGTH(secondset[0]));
		memcpy(Csort[i], Comp[i], sizeof(int)*LENGTH(secondset[0]));
		sort(Comp[i], Comp[i] + LENGTH(secondset[0]), compare);
	}
	int** CsortSeq = new int*[LENGTH(secondset)];
	for (int i = 0; i < LENGTH(secondset); i++)
	{
		CsortSeq[i] = new int[LENGTH(secondset[0])];
		memset(CsortSeq[i], 0, sizeof(int)*LENGTH(secondset[0]));
		//memcpy(CsortSeq[i], Comp[i], sizeof(int)*LENGTH(secondset[0]));
	}
	for (int i = 0; i < LENGTH(secondset); i++)
	{
		int tmpcol;
		for (tmpcol = 0; tmpcol < LENGTH(secondset[0]); tmpcol++)
		{
			for (int j = 0; j < LENGTH(secondset[0]); j++)
			{
				if (Csort[i][j] == Comp[i][tmpcol])
				{
					CsortSeq[i][tmpcol] = secondset[i][j];
					break;
				}
			}
		}
	}
	return CsortSeq;
}

//设置threshold
static double setThreshold(int ideal, int actual[]);

//资源是否被占用
static int isAvailable(MonitorPtr m, int obseq,int* _sj,int resid)
{
	int second = ABLE;
	int first = ABLE;
	for (int i = aj[obseq]; i <= bj[obseq]; i++)
	{
		if ((m->slot)[resid][i] == DISABLED) {
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
			if(m->slot[resid][j] == ABLE)
			cnt++;
		}
		if (cnt == pj[obseq]) {
			*_sj = j;
				return ABLE;
		}
	}
	return DISABLED;
}

//在给出的一组解中找出最好的解
int** findBestSolution(int** resmat,int* index,int* scnt)
{
	int no = LENGTH(resmat);
	int* numset = new int[no];
	for (int i = 0; i < no; i++) {
		numset[i] = getTracedNum(resmat[i]);
	}
	int max = *max_element(numset, numset + no);
	//int scnt = 0;
	*scnt = 0;
	index = new int[no];
	for (int i = 0; i < no; i++) {
		if (numset[i] == max) {
			scnt++;
			*index = i;
			index += 1;
		}
	}
	int** finalRes = new int*[*scnt];
	for (int j = 0; j < *scnt; j++)
	{
		finalRes[j] = new int[max];
		memcpy(finalRes[j], resmat[index[j]], sizeof(int)*max);
	}
	return finalRes;
}

static int getTracedNum(int* arr)
{
	int cnt = 0;
	for (int i = 0; i < LENGTH(arr); i++)
	{
		if (arr[i] == NOOBJECT)
			break;
		cnt++;
	}
	return cnt;
}

//计算探测中心位置
int** calculateCenterPos(int** targetList, int v[], int aj[], int bj[],int startP[][3],int** sj)
{

#if LINE_MODEL 

#endif
#if CIRCLE_MODEL 

#endif
}

//读Excel文件接口
void readExcel(char* filename)
{
	char* data[8];
	ifstream iFile(filename);
	if (!iFile.is_open())
	{
		errorReport("Can't open .csv file.");
		exit(-1);
	}
	char buffer[256];
	int i = 0,j = 0;
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
			//obarr[i - 1].curcor[1] = obarr[i - 1].curcor[2] = 0;
#endif	
		}
		j = 0;
		//cout << "Line " << i << ": " << buffer << endl;
		i++;
	}
}

//写Excel工作表接口
void writeExcel(int** res, int id) 
{

}

//将结果打印到终端
void printResult(int** res, int id)
{

}

//改写comp函数使sort从大到小排序
bool compare(const double &a, const double &b)
{
	return a>b;
}

//思考一个问题：每一次重置的时候，发生了什么？是不是每一次都更新了探测器的状态？
//逻辑全部写好之后再写到.h里面去
