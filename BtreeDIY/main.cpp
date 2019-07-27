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
//��ѯ��������
typedef map<colnum, pair<BtreeEntry, BtreeEntry>> getParams;

//��Ŷ����Զ��̵߳ĵ����̲߳�ѯ���
struct ResultSet
{
	//���ڱ�ʶ״̬��������ݼ��Ѿ�׼������Ϊtrue
	//Ҳ����sortByResultSetVec��������hash��������
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
	//��ʼ��hash��
	hashtable = (Listnode*)malloc(tablesize * sizeof(Listnode));
	for (sedKey i = 0; i < tablesize; i++){
		hashtable[i].nxtnodeptr = NULL;
		hashtable[i].killed = false;
		hashtable[i].pkeyhash = 0;
	}
}
void HashIntersector::finalResult(resultset & finalresult)
{
	/*����hash���õ����ս��finalresult*/
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
/*����hash��basecol�ǽ���hash��ĵ�һ�����Ե��������Ƕ��߳̽�����н������С���Ǹ������ı��*/
void HashIntersector::constructHashTable(mapResultSet&resultSet) {
	/*������*/
	vector<pKeyHash>::iterator it = resultSet[basecol]._result.begin();

	while (it!= resultSet[basecol]._result.end())
	{
		/*hash����H(key)=key MOD tablesize*/
		sedKey index = (*it) % tablesize;

		/*ָ����*/
		Listnode *ptr = hashtable[index].nxtnodeptr;
		/*ָ��ǰһ�����*/
		Listnode *qptr = hashtable+((*it) % tablesize);
		while (ptr!=NULL){
			qptr = ptr;
			ptr = ptr->nxtnodeptr;
		}
		/*�����½��*/
		Listnode *newnode = (Listnode *)malloc(sizeof(Listnode));
		/*��ʼ���½��*/
		(newnode)->nxtnodeptr = NULL;
		(newnode)->killed = true;
		(newnode)->pkeyhash = *it;
		/*���½������б�*/
		qptr->nxtnodeptr = newnode;

		itemnumber++;//��������һ
		it++;//������ǰ��
	}
}
/*����������*/
void HashIntersector::clean()
{
	for (sedKey i = 0; i < tablesize; i++)
	{
		Listnode *ptr=hashtable[i].nxtnodeptr;//ָ����
		Listnode *qptr = hashtable + i;//ָ��ǰһ�����
		while (ptr!=NULL)
		{
			if (ptr->killed){
				qptr->nxtnodeptr = ptr->nxtnodeptr;
				ptr->nxtnodeptr=NULL;
				free(ptr);//�ͷŽ��ռ�
				ptr = qptr->nxtnodeptr;
				itemnumber--;//��������һ
			}
			else{
				ptr->killed = true;//�ָ�killλ��Ϊ�´�ɸѡ��׼��
				//ָ��������
				ptr = ptr->nxtnodeptr;
				qptr = qptr->nxtnodeptr;
			}
			
		}
	}
}
void HashIntersector::seqHash(mapResultSet & resultSet)
{
	/*Ŀ�����Եı���*/
	for (colnum i = 1; i < resultSet.size(); i++)
	{
		/*��Ӧ���Ե����ݼ��ĵ����������ڱ������ݼ�*/
		vector<pKeyHash>::iterator it=(resultSet)[hashseq[i]]._result.begin();
		while (it != (resultSet)[hashseq[i]]._result.end())
		{
			/*ptrָ��������hash�еĽ��*/
			Listnode *ptr = (hashtable[(*it) % tablesize]).nxtnodeptr;
			while (ptr!=NULL)
			{
				/*���ָ��Ľ��ļ�ֵ�͵�ǰ�����ļ�ֵ��ͬ*/
				if (ptr->pkeyhash==(*it))
				{
					ptr->killed = false;//�ӱ���
				}
				ptr = ptr->nxtnodeptr;//ָ����һ�����
			}
			it++;
		}
		clean();//�����hash����killedΪtrue�Ľ��
	}
}
void HashIntersector::sortByResultSetVec(colnum * hashseq, mapResultSet & resultSet)
{
	//colnum min = NUM_TESTDATA;
	/*ֱ�Ӳ�������*/
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
		(resultSet)[basecol].state = false;//�ź���ı�־λfalse�����ٷ��ʡ�
		min = NUM_TESTDATA;
	}
}
HashIntersector::HashIntersector() {}
HashIntersector::~HashIntersector(){
	/*�ͷ�����Ŀռ䣬������ϣ��ͽ��*/
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
/*����hash�ṹ*/
void HASH_Intersector(mapResultSet& mapresultset, resultset&finalresult) {

	cout << "	*���ڽ���Hash�ṹ...";
	/*��ʼ��һ��HashIntersector��*/
	HashIntersector hashIntersector(mapresultset);
	/*����hash table*/
	hashIntersector.constructHashTable(mapresultset);
	/*���ݶ��̵߳õ������������������ν���hash��ײ����������HashIntersector���hashseq�У�*/
	hashIntersector.seqHash(mapresultset);
	/*���������յ�hashtable���õ����Ľ��*/
	hashIntersector.finalResult(finalresult);
	cout << "����ɡ�" << endl;
}
/*����ָ�����Ե�B+��*/
void initIndexBtree(
	std::map<pKeyHash, record>&srcdata
	, IndexBtree&btree
	, colnum col) {
	/*������*/
	std::map<pKeyHash, record>::iterator iter;
	iter = srcdata.begin();
	while (iter != srcdata.end())
	{
		BtreeEntry BE(iter->second.at(col), 100, iter->first);
		btree.insert(BE);
		iter++;
	}
}
/*�����������Ե�B+��*/
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
/*���ϸ�˹�ֲ��������*/
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
	param dtestdata:���ݼ�
	param size:��Ҫ���ɵ����ݼ���С
	param colnum:���ݼ���������
	param randparam:vector���ͽṹ��vector��СΪcolnum��С��ÿ�����������������Ե������������
*/
void dataGenerator(
	  std::map<pKeyHash, record>&dtestdata
	, uint64_t size
	, colnum colnum
	, uint16_t *col
	, std::vector<std::pair<double,double>> &randparam)
{
	if (colnum > randparam.size()) {
		std::cout<<"ERROR: �����������������"<<std::endl;
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
/*writeDataToFile����������д���ļ�*/
void writeDataToFile(
	  ofstream &oFile
	, const string filename
	, resultset &rs
	, std::map<pKeyHash, record>& testdata
	, colnum *col)
{
	cout << "	*���ڽ��������д���ļ�"<< filename <<"...";
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
	cout << "����ɡ�" << endl;
}
void writeDataToFile(
	  ofstream &oFile
	, const string filename
	, map<pKeyHash, record>&rs
	, std::map<pKeyHash, record>& testdata
	, colnum *col)
{
	cout << "	*���ڽ��������д���ļ�" << filename << "...";
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
	cout << "����ɡ�" << endl;
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
/*���ڶ��̲߳�ѯ���������Ľ����ÿ���߳�����һ����������Ҫ�����ǵ�����*/
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
/*���ڵ��Ʒ�ʽ�õ��������Ҫ�����ǵ�����*/
void recur_find(
	map<pKeyHash, record> &result
	, map<pKeyHash, record> &srcSet
	, pair<BtreeEntry, BtreeEntry> condition
	, IndexBtree *IndexBtree) {

	IndexBtree::iterator it_range = IndexBtree->lower_bound(condition.first);
	while (it_range != IndexBtree->upper_bound(condition.second))
	{
		result[it_range->pKeyHash]= srcSet[it_range->pKeyHash];//�õ���Ӧ��¼
		it_range++;
	}
}
/*�������߳�*/
void gettingSets(
	mapResultSet& results
	, colnum _colnum
	, colnum col[]
	, getParams& params
	, map<colnum, IndexBtree*>&IndexBtrees) {
	cout << "	*���ڽ��ж��߳�������ѯ...";
	bool flag=false;
	for (colnum i = 0; i < _colnum; i++)
	{
		thread t(hash_find, &results, col[i], params[col[i]],IndexBtrees[col[i]]);//col[i]�����Ա�ţ�params[col[i]]�ǲ�ѯ����pair<BtreeEntry, BtreeEntry>
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
	cout << "����ɡ�" << endl;
}
/*���������ǼǱ���������������*/
void indexDestory(map<colnum, IndexBtree*>&IndexBtrees, colnum *colindex, colnum colsize) {
	for (colnum i = 0; i < colsize; i++)
	{
		IndexBtrees[colindex[i]]->clear();
	}
}
void indexDestory(IndexBtree &IndexBtree) {
		IndexBtree.clear();
}
/*���Ʒ�ʽ*/
void Recur_Intersector(
	map<pKeyHash, record>& tempdata
	, map<pKeyHash, record> &testdata
	, colnum _colnum
	, colnum col[]
	, getParams& params
	, map<colnum, IndexBtree*>&IndexBtrees){

	cout << "	*���ڽ��ж��������Ʋ�ѯ...";
	tempdata.clear();
	IndexBtree tempIndexBtree(NUM_TESTATTRIBUTE);

	for (colnum i = 0; i < _colnum; i++)
	{
		/*�������*/
		tempIndexBtree.clear();
		if (i==0)
		{	/*��һ������ʹ�������ǼǱ��е������������ظ���������*/
			recur_find(tempdata, testdata, params[col[i]], IndexBtrees[col[i]]);
		}
		else
		{	/*���ݽ������������*/
			initIndexBtree(tempdata, tempIndexBtree, i);
			tempdata.clear();//��������
			//�����ϸ������������������ѯ���µĽ����
			recur_find(tempdata, testdata, params[col[i]], &tempIndexBtree);
		}
	}
	cout << "����ɡ�" << endl;
}
/*ȫɨ��*/
void fullScan(
	  map<pKeyHash, record> &outdata
	, map<pKeyHash, record> &indata
	, colnum _colnum
	, colnum col[]
	, getParams& params
) {
	cout << "	*���ڽ���ȫ��ɨ���ѯ...";
	map<pKeyHash, record>::iterator itera = indata.begin();

	while (itera != indata.end())
	{
		bool flag = true;
		/*̽��ÿ�������Ƿ���������*/
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
	cout << "����ɡ�" << endl;
}

int main() 
{
	/*
		���Գ���
		2019-3-15 ����
		���ֶ�ά���Բ�ѯ����
	*/
	//��������,ͬʱҲ��B+��������Ŀ��Ҳ�Ƿ���1�Ĳ����߳�����
	//uint16_t colnum = 0;
	setlocale(LC_ALL, "chs");
	ofstream oFile;
	cout << "�������ݼ���С��";
	cin >> NUM_TESTDATA;//�������ݼ���С
	cout <<endl;

	cout << "��������������";
	cin >> NUM_TESTATTRIBUTE;//����������
	cout << endl;
	
	/*		���������ݼ�	*/
	map<pKeyHash, record> testdata;
	/*	�������Ե���̬�ֲ�����	*/
	vector<pair<double, double>> randparam(NUM_MAX_TESTATTRIBUTE);	
	/*	��ʼ����̬�ֲ�����	*/
	cout << "���ڳ�ʼ����̬�ֲ�����...";
	randparam.at(0) = pair<double, double>(330, 100);/*��ֵ330������100����ͬ��*/
	randparam.at(1) = pair<double, double>(240, 100);
	randparam.at(2) = pair<double, double>(300, 100);
	randparam.at(3) = pair<double, double>(180, 100);
	randparam.at(4) = pair<double, double>(360, 100);
	randparam.at(5) = pair<double, double>(210, 100);
	randparam.at(6) = pair<double, double>(150, 100);
	randparam.at(7) = pair<double, double>(270, 100);
	/*�������ݣ�NUM_TESTDATA��������NUM_TESTATTRIBUTE����������*/
	cout << "����ɡ�" << endl;
	/*���ڴ洢��Ҫ�������������Ա��*/
	colnum *colindex=(colnum*)malloc(NUM_TESTATTRIBUTE*sizeof(colnum));

	cout << "���뽨��������"<<NUM_TESTATTRIBUTE<<"�����Ե��к�[0,7]:";
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++){
		cin>> colindex[i];//���ܿ���̨��������Ա��
	}
	cout << "�������ɲ�������...";
	/*���ɲ������ݵĺ���*/
	dataGenerator(testdata, NUM_TESTDATA, NUM_TESTATTRIBUTE, colindex, randparam);
	cout << "����ɡ�" << endl;
	
	/*	�����������б�ͨ���˿����ҵ����н����������ṹ	*/
	cout << "���ڳ�ʼ�������ǼǱ�...";
	std::map<colnum, IndexBtree*> indexBtrees;
	/*	��ʼ�������б�	*/
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++){
		IndexBtree *btree = new IndexBtree(colindex[i]);//�����ṹ���б��
		indexBtrees[colindex[i]] = btree;
	}
	cout << "����ɡ�" << endl;

	/*Ϊ������������������ B+ ������*/
	cout << "���ڳ�ʼ������...";
	initIndexBtree(testdata, indexBtrees, colindex);
	cout << "����ɡ�" << endl;
	
	mapResultSet hashResultSet;
	getParams params;

	cout << "�������ò�ѯ����..."<<endl;
	for (colnum i = 0; i < NUM_TESTATTRIBUTE; i++) {
		sedKey low, high;
		cout << "����" << colindex[i] << "�����ԵĲ�ѯ��Χ��" << endl;
		cin >> low;
		cin >> high;
		params[colindex[i]] = pair<BtreeEntry, BtreeEntry >(BtreeEntry(low, 100, 0), BtreeEntry(high, 100, 0));
	}
	cout << "����ɡ�" << endl;

	cout << "*****************���߳�Hash����*******************" << endl;
	DWORD StartTime = ::GetTickCount();
	gettingSets(hashResultSet, NUM_TESTATTRIBUTE, colindex, params, indexBtrees);
	vector<pKeyHash> hashrs;//��Ž��
	HASH_Intersector(hashResultSet,hashrs);
	DWORD EndTime = ::GetTickCount();
	cout << "	*��ʱ" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile,"result_hash.csv", hashrs, testdata,colindex);
	

	cout << "*****************���������Ʒ���*******************" << endl;
	StartTime = ::GetTickCount();
	map<pKeyHash, record> recurrs;//��Ž��
	Recur_Intersector(recurrs,testdata ,NUM_TESTATTRIBUTE, colindex, params, indexBtrees);
	EndTime = ::GetTickCount();
	cout << "	*��ʱ" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile, "result_recur.csv", recurrs, testdata, colindex);

	/***************************************/
	cout << "*****************ȫ��˳��鷽��*******************" << endl;

	StartTime = ::GetTickCount();
	map<pKeyHash, record> fullscrc;//��Ž��
	fullScan(fullscrc, testdata, NUM_TESTATTRIBUTE, colindex, params);
	EndTime = ::GetTickCount();
	cout << "	*��ʱ" << EndTime - StartTime << "ms" << endl;
	writeDataToFile(oFile, "result_fulls.csv", fullscrc, testdata, colindex);

	cout << "**************************************************" << endl;
	/*����������д���ļ�*/
	cout << "���ڽ���������д���ļ�...";
	writeDataToFile(oFile,testdata, colindex);
	cout << "����ɡ�" << endl;

	cout << "������������...";
	indexDestory(indexBtrees, colindex, NUM_TESTATTRIBUTE);
	free(colindex);
	cout << "����ɡ�" << endl;
	
	getchar();
	getchar();
}
