#include <map>
#include <cmath>
#include <vector>
#include <limits>
# include<thread>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <windows.h>
#include "btree.h"

//#define NUM_TESTATTRIBUTE 4
#define NUM_MAX_TESTATTRIBUTE 8
//#define NUM_TESTDATA 100000

using namespace std;

typedef uint64_t pKeyHash;
typedef uint16_t colnum;
typedef vector<sedKey> record;

typedef vector<pKeyHash> resultset;
//查询条件参数
typedef map<colnum, pair<BtreeEntry, BtreeEntry>> getParams;

//存放多属性多线程的单个线程查询结果
struct ResultSet
{
	//用于标识状态，如果数据集已经准备好则为true
	//也用于sortByResultSetVec函数进行hash序列排序
	bool state=false;		
	resultset _result;	
	void reset() {
		state = false;
		_result.clear();
	}
};

typedef map<colnum, ResultSet> mapResultSet;
sedKey NUM_TESTDATA = 0;
colnum NUM_TESTATTRIBUTE = 0;
class HashIntersector
{
public:
	HashIntersector();
	HashIntersector(mapResultSet&resultSet);
	void finalResult(resultset& finalresult);
	void constructHashTable(mapResultSet&resultSet);
	void seqHash(mapResultSet &resultSet);
	~HashIntersector();

private:

	struct Listnode{
		Listnode* nxtnodeptr=NULL;
		bool killed=false;
		pKeyHash pkeyhash=0;
	};

	sedKey itemnumber = 0;
	colnum basecol=0;
	size_t tablesize=0;
	colnum* hashseq = NULL;
	Listnode* hashtable=NULL;

	void clean();
	void sortByResultSetVec(colnum *hashseq, mapResultSet &resultSet);
};

HashIntersector::HashIntersector(mapResultSet&resultSet)
{
	hashseq = (colnum*)malloc(resultSet.size() * sizeof(colnum));
	
	sortByResultSetVec(hashseq, resultSet);
	if (hashseq)
	{
		basecol = hashseq[0];
	}

	tablesize = resultSet[basecol]._result.size();
	//初始化hash表
	hashtable = (Listnode*)malloc(tablesize * sizeof(Listnode));
	for (sedKey i = 0; i < tablesize; i++){
		hashtable[i].nxtnodeptr = NULL;
		hashtable[i].killed = false;
		hashtable[i].pkeyhash = 0;
	}
}
void HashIntersector::finalResult(resultset & finalresult)
{
	/*遍历hash表，得到最终结果finalresult*/
	for (sedKey i = 0; i < tablesize; i++)
	{
		Listnode *ptr = hashtable[i].nxtnodeptr;
		Listnode *qptr = hashtable + i;
		while (ptr != NULL)
		{
			finalresult.push_back(ptr->pkeyhash);
			ptr = ptr->nxtnodeptr;
			qptr = qptr->nxtnodeptr;
		}
	}
}
/*建立hash表，basecol是建立hash表的第一个属性的索引，是多线程结果集中结果数最小的那个索引的编号*/
void HashIntersector::constructHashTable(mapResultSet&resultSet) {
	/*迭代器*/
	vector<pKeyHash>::iterator it = resultSet[basecol]._result.begin();

	while (it!= resultSet[basecol]._result.end())
	{
		/*hash函数H(key)=key MOD tablesize*/
		sedKey index = (*it) % tablesize;

		/*指向结点*/
		Listnode *ptr = hashtable[index].nxtnodeptr;
		/*指向前一个结点*/
		Listnode *qptr = hashtable+((*it) % tablesize);
		while (ptr!=NULL){
			qptr = ptr;
			ptr = ptr->nxtnodeptr;
		}
		/*申请新结点*/
		Listnode *newnode = (Listnode *)malloc(sizeof(Listnode));
		/*初始化新结点*/
		(newnode)->nxtnodeptr = NULL;
		(newnode)->killed = true;
		(newnode)->pkeyhash = *it;
		/*将新结点接入列表*/
		qptr->nxtnodeptr = newnode;

		itemnumber++;//计数器加一
		it++;//迭代器前进
	}
}
/*清除废弃结点*/
void HashIntersector::clean()
{
	for (sedKey i = 0; i < tablesize; i++)
	{
		Listnode *ptr=hashtable[i].nxtnodeptr;//指向结点
		Listnode *qptr = hashtable + i;//指向前一个结点
		while (ptr!=NULL)
		{
			if (ptr->killed){
				qptr->nxtnodeptr = ptr->nxtnodeptr;
				ptr->nxtnodeptr=NULL;
				free(ptr);//释放结点空间
				ptr = qptr->nxtnodeptr;
				itemnumber--;//计数器减一
			}
			else{
				ptr->killed = true;//恢复kill位，为下次筛选做准备
				//指向下组结点
				ptr = ptr->nxtnodeptr;
				qptr = qptr->nxtnodeptr;
			}
			
		}
	}
}
void HashIntersector::seqHash(mapResultSet & resultSet)
{
	/*目标属性的遍历*/
	for (colnum i = 1; i < resultSet.size(); i++)
	{
		/*对应属性的数据集的迭代器，用于遍历数据集*/
		vector<pKeyHash>::iterator it=(resultSet)[hashseq[i]]._result.begin();
		while (it != (resultSet)[hashseq[i]]._result.end())
		{
			/*ptr指向拉链法hash中的结点*/
			Listnode *ptr = (hashtable[(*it) % tablesize]).nxtnodeptr;
			while (ptr!=NULL)
			{
				/*如果指向的结点的键值和当前遍历的键值相同*/
				if (ptr->pkeyhash==(*it))
				{
					ptr->killed = false;//加保护
				}
				ptr = ptr->nxtnodeptr;//指向下一个结点
			}
			it++;
		}
		clean();//清理掉hash表中killed为true的结点
	}
}
void HashIntersector::sortByResultSetVec(colnum * hashseq, mapResultSet & resultSet)
{
	//colnum min = NUM_TESTDATA;
	/*直接插入排序*/
	mapResultSet::iterator it = (resultSet).begin();

	sedKey min = NUM_TESTDATA;
	for (colnum i = 0; i < (resultSet).size(); i++){
		mapResultSet::iterator tempit = (resultSet).begin();
		while (tempit!= (resultSet).end())
		{
			if ((tempit->second.state==true)&&tempit->second._result.size() <= min) {
				basecol = tempit->first;
				min = tempit->second._result.size();
			}
			tempit++;
		}
		hashseq[i] = basecol;
		(resultSet)[basecol].state = false;//排好序的标志位false，不再访问。
		min = NUM_TESTDATA;
	}
}
HashIntersector::HashIntersector() {}
HashIntersector::~HashIntersector(){
	/*释放申请的空间，包括哈希表和结点*/
	for (sedKey i = 0; i < tablesize; i++){
		if (hashtable[i].nxtnodeptr != NULL) {
			Listnode* ptr = hashtable[i].nxtnodeptr;
			while (ptr){
				Listnode* q = ptr;
				ptr = ptr->nxtnodeptr;
				free(q);
			}
		}
	}
	if (hashtable != NULL){
		free(hashtable);
	}
	if (hashseq != NULL) {
		free(hashseq);
	}
}
/*建立hash结构*/
void HASH_Intersector(mapResultSet& mapresultset, resultset&finalresult) {

	cout << "	*正在建立Hash结构...";
	/*初始化一个HashIntersector类*/
	HashIntersector hashIntersector(mapresultset);
	/*建立hash table*/
	hashIntersector.constructHashTable(mapresultset);
	/*根据多线程得到的数据量排序结果依次进行hash碰撞（排序结果在HashIntersector类的hashseq中）*/
	hashIntersector.seqHash(mapresultset);
	/*遍历做最终的hashtable，得到最后的结果*/
	hashIntersector.finalResult(finalresult);
	cout << "【完成】" << endl;
}
/*建立指定属性的B+树*/
void initIndexBtree(
	std::map<pKeyHash, record>&srcdata
	, IndexBtree&btree
	, colnum col) {
	/*迭代器*/
	std::map<pKeyHash, record>::iterator iter;
	iter = srcdata.begin();
	while (iter != srcdata.end())
	{
		BtreeEntry BE(iter->second.at(col), 100, iter->first);
		btree.insert(BE);
		iter++;
	}
}
/*建立所有属性的B+树*/
void initIndexBtree(
	std::map<pKeyHash, record>&srcdata
	, std::map<colnum, IndexBtree*>&indexBtrees
	, colnum col[]) {

	std::map<colnum, IndexBtree*>::iterator it;
	it = indexBtrees.begin();

	for (colnum i = 0; i < indexBtrees.size(); i++)
	{
		std::map<pKeyHash, record>::iterator iter;
		iter = srcdata.begin();
		while (iter != srcdata.end())
		{
			//BtreeEntry BE(iter->second.at(it->first), 100, iter->first);
			BtreeEntry BE(iter->second.at(i), 100, iter->first);
			indexBtrees[col[i]]->insert(BE);
			iter++;
		}
	}


	//while (it != indexBtrees.end())
	//{
	//	std::map<pKeyHash, record>::iterator iter;
	//	iter = srcdata.begin();
	//	while (iter != srcdata.end())
	//	{

	//		//if (it->first < iter->second.size()) {
	//			//BtreeEntry *BE=new BtreeEntry(iter->second.at(it->first), 100, iter->first);
	//			BtreeEntry BE(iter->second.at(it->first), 100, iter->first);
	//			it->second->insert(BE);
	//		//}
	//		//else
	//		/*{
	//			std::cout << "the column number is out of the attributes' number! " << std::endl;
	//			break;
	//		}*/
	//		iter++;
	//	}
	//	it++;
	//}
}
/*符合高斯分布的随机数*/
double generateGaussianNoise(double mu, double sigma)
{
	const double epsilon = 2.2250738585072014e-308;
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	} while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}
/*
	param dtestdata:数据集
	param size:想要生成的数据集大小
	param colnum:数据集的属性数
	param randparam:vector类型结构，vector大小为colnum大小，每个分量代表所在属性的随机函数参数
*/
void dataGenerator(
	  std::map<pKeyHash, record>&dtestdata
	, uint64_t size
	, colnum colnum
	, uint16_t *col
	, std::vector<std::pair<double,double>> &randparam)
{
	if (colnum > randparam.size()) {
		std::cout<<"ERROR: 超出最大属性数限制"<<std::endl;
		return;
	}
	dtestdata.clear();
	for (uint64_t i = 0; i < size; i++)
	{
		record records(colnum);
		for (uint16_t j = 0; j < colnum; j++)
		{
			records.at(j) = sedKey(abs(generateGaussianNoise(randparam.at(col[j]).first, randparam.at(col[j]).second)));
		}
		dtestdata[i] = records;
	}

}
/*writeDataToFile函数将数据写入文件*/
void writeDataToFile(
	  ofstream &oFile
	, const string filename
	, resultset &rs
	, std::map<pKeyHash, record>& testdata
	, colnum *col)
{
	cout << "	*正在将结果数据写入文件"<< filename <<"...";
	oFile.open(filename, ios::out | ios::trunc);

	vector<pKeyHash>::iterator itera = rs.begin();
	oFile << "pkeyHash";
	for (int i = 0; i < testdata.at(i).size(); i++)
	{
		oFile << ",col" << col[i];
	}
	oFile << std::endl;
	while (itera != rs.end())
	{
		oFile << *itera;
		for (int i = 0; i < testdata[*itera].size(); i++)
		{
			oFile << "," << (testdata[*itera]).at(i);
		}
		oFile << std::endl;
		itera++;
	}
	oFile.close();
	cout << "【完成】" << endl;
}
void writeDataToFile(
	  ofstream &oFile
	, const string filename
	, map<pKeyHash, record>&rs
	, std::map<pKeyHash, record>& testdata
	, colnum *col)
{
	cout << "	*正在将结果数据写入文件" << filename << "...";
	oFile.open(filename, ios::out | ios::trunc);

	map<pKeyHash, record>::iterator itera = rs.begin();
	oFile << "pkeyHash";
	for (int i = 0; i < testdata.at(i).size(); i++)
	{
		oFile << ",col" << col[i];
	}
	oFile << std::endl;
	while (itera != rs.end())
	{
		oFile << itera->first;
		for (int i = 0; i < itera->second.size(); i++)
		{
			oFile << "," << itera->second.at(i);
		}
		oFile << std::endl;
		itera++;
	}
	oFile.close();
	cout << "【完成】" << endl;
}
void writeDataToFile(
	  ofstream &oFile
	, std::map<pKeyHash, record>& testdata
	, colnum *col)
{
	oFile.open("testdata.csv", ios::out | ios::trunc);

	map<pKeyHash, record>::iterator itera = testdata.begin();
	oFile << "pkeyHash";
	for (int i = 0; i < testdata.at(i).size(); i++)
	{
		oFile << ",col" << col[i];
	}
	oFile << std::endl;
	while (itera!= testdata.end())
	{
		oFile << itera->first ;
		for (int i=0;i< testdata.at(i).size();i++)
		{
			oFile <<"," << itera->second.at(i);
		}
		oFile << std::endl;
		itera++;
	}
	oFile.close();
}
/*用于多线程查询满足条件的结果，每个线程运行一个函数，主要工具是迭代器*/
void hash_find(
	  mapResultSet* results
	, colnum col
	, pair<BtreeEntry, BtreeEntry> condition
	, IndexBtree *IndexBtree) {

	IndexBtree::iterator it_range = IndexBtree->lower_bound(condition.first);
	while (it_range != IndexBtree->upper_bound(condition.second))
	{
		(*results)[col]._result.push_back(it_range->pKeyHash);
		it_range++;
	}
	(*results)[col].state = true;
}
/*用于递推方式得到结果，主要工具是迭代器*/
void recur_find(
	map<pKeyHash, record> &result
	, map<pKeyHash, record> &srcSet
	, pair<BtreeEntry, BtreeEntry> condition
	, IndexBtree *IndexBtree) {

	IndexBtree::iterator it_range = IndexBtree->lower_bound(condition.first);
	while (it_range != IndexBtree->upper_bound(condition.second))
	{
		result[it_range->pKeyHash]= srcSet[it_range->pKeyHash];//得到对应记录
		it_range++;
	}
}
/*启动多线程*/
void gettingSets(
	mapResultSet& results
	, colnum _colnum
	, colnum col[]
	, getParams& params
	, map<colnum, IndexBtree*>&IndexBtrees) {
	cout << "	*正在进行多线程条件查询...";
	bool flag=false;
	for (colnum i = 0; i < _colnum; i++)
	{
		thread t(hash_find, &results, col[i], params[col[i]],IndexBtrees[col[i]]);//col[i]是属性编号，params[col[i]]是查询条件pair<BtreeEntry, BtreeEntry>
		t.detach();
	}

	while (!flag)
	{
		bool tempflag = true;
		for (colnum i = 0; i < _colnum; i++)
		{
			tempflag &= results[col[i]].state;
		}
		flag = tempflag;
	}
	cout << "【完成】" << endl;
}
/*销毁索引登记表中所有索引索引*/
void indexDestory(map<colnum, IndexBtree*>&IndexBtrees, colnum *colindex, colnum colsize) {
	for (colnum i = 0; i < colsize; i++)
	{
		IndexBtrees[colindex[i]]->clear();
	}
}
void indexDestory(IndexBtree &IndexBtree) {
		IndexBtree.clear();
}
/*递推方式*/
void Recur_Intersector(
	map<pKeyHash, record>& tempdata
	, map<pKeyHash, record> &testdata
	, colnum _colnum
	, colnum col[]
	, getParams& params
	, map<colnum, IndexBtree*>&IndexBtrees){

	cout << "	*正在进行多索引递推查询...";
	tempdata.clear();
	IndexBtree tempIndexBtree(NUM_TESTATTRIBUTE);

	for (colnum i = 0; i < _colnum; i++)
	{
		/*清除索引*/
		tempIndexBtree.clear();
		if (i==0)
		{	/*第一个属性使用索引登记表中的索引，避免重复创建索引*/
			recur_find(tempdata, testdata, params[col[i]], IndexBtrees[col[i]]);
		}
		else
		{	/*根据结果集建立索引*/
			initIndexBtree(tempdata, tempIndexBtree, i);
			tempdata.clear();//清除结果集
			//根据上个结果集构建的索引查询出新的结果集
			recur_find(tempdata, testdata, params[col[i]], &tempIndexBtree);
		}
	}
	cout << "【完成】" << endl;
}
/*全扫描*/
void fullScan(
	  map<pKeyHash, record> &outdata
	, map<pKeyHash, record> &indata
	, colnum _colnum
	, colnum col[]
	, getParams& params
) {
	cout << "	*正在进行全表扫描查询...";
	map<pKeyHash, record>::iterator itera = indata.begin();

	while (itera != indata.end())
	{
		bool flag = true;
		/*探测每个属性是否满足条件*/
		for (colnum i = 0; i < _colnum; i++){
			if ((itera->second.at(i) >= params[col[i]].first.sKey) && (itera->second.at(i) < params[col[i]].second.sKey))
				flag &= true;
			else
				flag &= false;
		}
		if (flag){
			outdata[itera->first] = itera->second;
		}
		itera++;
	}
	cout << "【完成】" << endl;
}

int main() 
{
	/*
		测试程序
		2019-3-15 安迪
		两种多维属性查询方法
	*/
	//属性列数,同时也是B+树索引数目，也是方法1的查找线程数；
	//uint16_t colnum = 0;
	setlocale(LC_ALL, "chs");
	ofstream oFile;
	cout << "设置数据集大小：";
	cin >> NUM_TESTDATA;//输入数据集大小
	cout <<endl;

	cout << "设置属性列数：";
	cin >> NUM_TESTATTRIBUTE;//输入属性数
	cout << endl;
	
	/*		测试用数据集	*/
	map<pKeyHash, record> testdata;
	/*	各个属性的正态分布参数	*/
	vector<pair<double, double>> randparam(NUM_MAX_TESTATTRIBUTE);	
	/*	初始化正态分布参数	*/
	cout << "正在初始化正态分布参数...";
	randparam.at(0) = pair<double, double>(330, 100);/*均值330，方差100（下同）*/
	randparam.at(1) = pair<double, double>(240, 100);
	randparam.at(2) = pair<double, double>(300, 100);
	randparam.at(3) = pair<double, double>(180, 100);
	randparam.at(4) = pair<double, double>(360, 100);
	randparam.at(5) = pair<double, double>(210, 100);
	randparam.at(6) = pair<double, double>(150, 100);
	randparam.at(7) = pair<double, double>(270, 100);
	/*生成数据，NUM_TESTDATA数据量；NUM_TESTATTRIBUTE测试属性量*/
	cout << "【完成】" << endl;
	/*用于存储将要建立索引的属性编号*/
	colnum *colindex=(colnum*)malloc(NUM_TESTATTRIBUTE*sizeof(colnum));

	cout << "输入建立索引的"<<NUM_TESTATTRIBUTE<<"个属性的列号[0,7]:";
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++){
		cin>> colindex[i];//接受控制台输入的属性编号
	}
	cout << "正在生成测试数据...";
	/*生成测试数据的函数*/
	dataGenerator(testdata, NUM_TESTDATA, NUM_TESTATTRIBUTE, colindex, randparam);
	cout << "【完成】" << endl;
	
	/*	属性索引的列表，通过此可以找到所有建立的索引结构	*/
	cout << "正在初始化索引登记表...";
	std::map<colnum, IndexBtree*> indexBtrees;
	/*	初始化索引列表	*/
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++){
		IndexBtree *btree = new IndexBtree(colindex[i]);//索引结构进行编号
		indexBtrees[colindex[i]] = btree;
	}
	cout << "【完成】" << endl;

	/*为数据所有属性列生成 B+ 树索引*/
	cout << "正在初始化索引...";
	initIndexBtree(testdata, indexBtrees, colindex);
	cout << "【完成】" << endl;
	
	mapResultSet hashResultSet;
	getParams params;

	cout << "正在设置查询条件..."<<endl;
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++) {
		sedKey low, high;
		cout << "设置" << colindex[i] << "号属性的查询范围：" << endl;
		cin >> low;
		cin >> high;
		params[colindex[i]] = pair<BtreeEntry, BtreeEntry >(BtreeEntry(low, 100, 0), BtreeEntry(high, 100, 0));
	}
	cout << "【完成】" << endl;

	cout << "*****************多线程Hash方法*******************" << endl;
	DWORD StartTime = ::GetTickCount();
	gettingSets(hashResultSet, NUM_TESTATTRIBUTE, colindex, params, indexBtrees);
	vector<pKeyHash> hashrs;//存放结果
	HASH_Intersector(hashResultSet,hashrs);
	DWORD EndTime = ::GetTickCount();
	cout << "	*耗时" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile,"result_hash.csv", hashrs, testdata,colindex);
	

	cout << "*****************多索引递推方法*******************" << endl;
	StartTime = ::GetTickCount();
	map<pKeyHash, record> recurrs;//存放结果
	Recur_Intersector(recurrs,testdata ,NUM_TESTATTRIBUTE, colindex, params, indexBtrees);
	EndTime = ::GetTickCount();
	cout << "	*耗时" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile, "result_recur.csv", recurrs, testdata, colindex);

	/***************************************/
	cout << "*****************全表顺序查方法*******************" << endl;

	StartTime = ::GetTickCount();
	map<pKeyHash, record> fullscrc;//存放结果
	fullScan(fullscrc, testdata, NUM_TESTATTRIBUTE, colindex, params);
	EndTime = ::GetTickCount();
	cout << "	*耗时" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile, "result_fulls.csv", fullscrc, testdata, colindex);

	cout << "**************************************************" << endl;
	/*将测试数据写入文件*/
	cout << "正在将测试数据写入文件...";
	writeDataToFile(oFile,testdata, colindex);
	cout << "【完成】" << endl;

	cout << "正在销毁索引...";
	indexDestory(indexBtrees, colindex, NUM_TESTATTRIBUTE);
	free(colindex);
	cout << "【完成】" << endl;
	
	getchar();
	getchar();
}
