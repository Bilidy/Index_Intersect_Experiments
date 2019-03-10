#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>

#define INVALID_NODEID (1)
#define ROOT_ID 100

typedef int int32_t;
typedef unsigned short int uint16_t;
typedef unsigned int uint32_t;
typedef uint64_t NodeId;
typedef uint32_t sedKey;

struct BtreeEntry
{
	sedKey sKey;
	uint16_t sKeyLength;
	uint64_t pKeyHash;

	BtreeEntry(sedKey _sKey, uint16_t _sKeyLength, uint64_t _pKeyHash)
		:sKey(_sKey)
		, sKeyLength(_sKeyLength)
		, pKeyHash(_pKeyHash) {}
	BtreeEntry(sedKey _sKey, uint64_t _pKeyHash)
		: sKey(_sKey)
		, sKeyLength(sizeof(_sKey))
		, pKeyHash(_pKeyHash) {}
	BtreeEntry()
		: sKey(NULL)
		, sKeyLength(0)
		, pKeyHash(0) {}
	BtreeEntry& operator = (const BtreeEntry& other) 
	{
		this->sKey = other.sKey;
		this->sKeyLength = other.sKeyLength;
		this->pKeyHash = other.pKeyHash;
		return *this;
	}
	bool operator == (const BtreeEntry& other)
	{
		return
			(sKeyLength == other.sKeyLength) &&
			(pKeyHash == other.pKeyHash)&&
			(sKey == other.sKey);
	}
	bool operator != (const BtreeEntry& other)
	{
		return !(operator == (other));
	}
};
class IndexBtree
{
public:
	static const uint16_t leafslotmax = 8;
	static const uint16_t innerslotmax = leafslotmax;
	static const uint16_t minleafslots = (leafslotmax / 2);
	static const uint16_t mininnerslots = (innerslotmax / 2);
	static const bool selfverify = false;
	static const bool useBinarySearch = true;
	struct tree_stats {
		uint64_t itemcount;// ���ٸ�items
		uint64_t leaves;// ���ٲ�
		uint64_t innernodes;// �ڽڵ���
		static const uint16_t leafslots = leafslotmax;// ÿ��Ҷ�ڵ��slot�� leafslotmax = 8
		static const uint16_t innerslots = innerslotmax;// ÿ���ڽڵ��slot�� innerslotmax = leafslotmax
		// minleafslots = 4;  leafslotmax = leafslots = 8
		// mininnerslots = 4; innerslotmax = innerslots = 8 
		inline tree_stats()
			: itemcount(0),
			leaves(0),
			innernodes(0)
		{}
		inline uint64_t nodes() const {// ������node��
			return innernodes + leaves;
		}
		inline double avgfill_leaves() const {// Return the average fill of leaves
			return static_cast<double>((itemcount) / (leaves * leafslots));
		}
	};
	

	IndexBtree();
	~IndexBtree();

//private:
	tree_stats  m_stats;
	uint64_t treeTableId;
	uint64_t nextNodeId;
	uint32_t numEntries;
	NodeId m_rootId;

	typedef
	struct Node {

		struct KeyInfo {
			sedKey sKey=NULL;
			//int32_t relOffset;
			uint16_t keyLength;
			uint64_t pkHash;
			/*uint32_t endRelOffset() {
				return relOffset + keyLength;
			}*/
			KeyInfo() :keyLength(0), pkHash(0) 
			{
			};
			void clear(){
				sKey = NULL;
				keyLength = 0;
				pkHash = 0;
			}
			KeyInfo& operator=(const KeyInfo&other)
			{
				this->sKey = other.sKey;
				this->pkHash = other.pkHash;
				this->keyLength = other.keyLength;
				return *this;
			}
		};
		uint32_t keysBeginOffset;
		uint16_t level;
		uint16_t slotuse;
		uint32_t keyStorageUsed;
		KeyInfo keys[IndexBtree::innerslotmax];
		/*Node(char *backingStorage, uint16_t level = 0)
			: keyBuffer(backingStorage)
			, keysBeginOffset(strlen(backingStorage))
			, level(level)
			, slotuse(0)
			, keyStorageUsed(0)
		{}*/
		Node(uint16_t level = 0)
			: level(level)
			, slotuse(0)
			, keyStorageUsed(0) {}
		virtual ~Node() {}
		inline BtreeEntry 
			getAt(uint16_t index) const
		{
			if (!(index < Node::slotuse)) 
			{
				printf("getAt() erorr:BtreeEntry not exist!\n exiting\n");
				exit(0);
			}
			return BtreeEntry(keys[index].sKey, keys[index].keyLength, keys[index].pkHash);
		};
		inline void 
			setAt(uint16_t index,const BtreeEntry entry)
		{
			if (!(index < IndexBtree::innerslotmax))
			{
				printf("setAt() erorr:BtreeEntry index out of bound!\n exiting\n");
				exit(0);
			};
			if (index < Node::slotuse)
			{
				keys[index].sKey = entry.sKey;
				keys[index].keyLength = entry.sKeyLength;
				keys[index].pkHash = entry.pKeyHash;
			}
			else
			{
				keys[index].sKey = entry.sKey;
				keys[index].keyLength = entry.sKeyLength;
				keys[index].pkHash = entry.pKeyHash;
				slotuse++;
			}
		}
		inline BtreeEntry 
			back() const
		{
			if (0 < Node::slotuse)
			{
				return getAt(uint16_t(slotuse - 1));
			}
		}
		inline void
			pop_back(uint16_t n = 1)
		{
			if (Node::slotuse >= n) 
			{
				printf("pop_back() erorr:The number of poped is greater than soltuse!\n exiting\n");
				exit(0);
			};

			// Do a virtual copy instead of just decrementing buffer size so
			// any other node's references to the old data will still be valid
			slotuse = uint16_t(slotuse - n);
		}
		inline bool 
			isLeaf() const{
			return (level == 0);
		}
		inline bool 
			isinnenode() const {
			return !isLeaf();
		}
		inline bool 
			isfull()const {
			return (Node::slotuse == innerslotmax);
		}
		inline bool 
			isfew() const{
			return (Node::slotuse <= mininnerslots);
		}
		inline bool 
			isunderflow() const{
			return (Node::slotuse < mininnerslots);
		}
		//protected:

		inline void 
			insertAtEntryOnly(uint16_t index, BtreeEntry entry){
			if ((index > Node::slotuse))
			{
				printf("insertAtEntryOnly() erorr:BtreeEntry not exist!\n exiting\n");
				exit(0);
			}
			if (index < Node::slotuse)//index=3 slotuse=4
			{
				for (uint32_t i = slotuse; i >= index + 1; i--)
				{
					keys[i] = keys[i - 1];
				}
				keys[index].sKey = entry.sKey;
				keys[index].pkHash = entry.pKeyHash;
				keys[index].keyLength = entry.sKeyLength;
				slotuse++;
			}
			else//index=Node::slotuse index=4 slotuse=4
			{
				keys[index].sKey = entry.sKey;
				keys[index].pkHash = entry.pKeyHash;
				keys[index].keyLength = entry.sKeyLength;
				slotuse++;
			}
			
		}
		inline void 
			eraseAtEntryOnly(uint16_t index){
			if (index > Node::slotuse)
			{
				printf("eraseAtEntryOnly() erorr:index out of bound!\n exiting\n");
				exit(0);
			}
			if (index < Node::slotuse)
			{
				for (uint16_t i = index; i < Node::slotuse; i++)
				{
					keys[i] = keys[i+1];
				}
				keys[Node::slotuse].clear();
				slotuse--;
			}
			else//index = Node::slotuse
			{
				keys[Node::slotuse].clear();
				slotuse--;
			}
			
		}
		inline void
			moveBackEntriesToFrontOf(Node *dest, uint16_t numEntries) {
			if (numEntries + dest->slotuse > IndexBtree::innerslotmax) 
			{
				printf("moveBackEntriesToFrontOf() erorr:dest BtreeEntry's slotuse will greater than innerslotmax!\n exiting\n");
				exit(0);
			}
			if (numEntries > slotuse)
			{
				printf("moveBackEntriesToFrontOf() erorr:the number of Entries is greater than src's slotuse!\n exiting\n");
				exit(0);
			}
			uint16_t splitPoint = uint16_t(slotuse - numEntries);
			for (uint16_t i = Node::slotuse - 1; i >= splitPoint; i--)
			{
				BtreeEntry BE(keys[i].sKey, keys[i].keyLength, keys[i].pkHash);
				dest->insertAtEntryOnly(0, BE);
				eraseAtEntryOnly(i);
			}
		}
		inline void
			moveFrontEntriesToBackOf(Node *dest, uint16_t numEntries) {
			if(numEntries + dest->slotuse> IndexBtree::innerslotmax){
				printf("moveFrontEntriesToBackOf() erorr:dest BtreeEntry's slotuse will greater than innerslotmax!\n exiting\n");
				exit(0);
			}
			if(numEntries > slotuse){
				printf("moveFrontEntriesToBackOf() erorr:the number of Entries is greater than src's slotuse!\n exiting\n");
				exit(0);
			}
			uint16_t splitPoint = uint16_t(numEntries - 1);
			for (uint16_t i = 0; i <= splitPoint; i++)
			{
				BtreeEntry BE(keys[0].sKey, keys[0].keyLength, keys[0].pkHash);
				dest->insertAtEntryOnly(dest->slotuse, BE);
				eraseAtEntryOnly(0);
			}

		}
	}Node, *Nodeptr;
	struct InnerNode : public Node 
	{
		NodeId child[IndexBtree::innerslotmax + 1];
		KeyInfo rightMostLeafKey;
		bool rightMostLeafKeyIsInfinite;

		InnerNode(uint16_t level)
			: child()
			, rightMostLeafKey()
			, rightMostLeafKeyIsInfinite(true){}

		void setRightMostLeafKey(BtreeEntry entry) {

			rightMostLeafKey.sKey = entry.sKey;
			rightMostLeafKey.keyLength = entry.sKeyLength;
			rightMostLeafKey.pkHash = entry.pKeyHash;

			rightMostLeafKeyIsInfinite = false;
		}

		BtreeEntry getRightMostLeafKey() {
			return { 
				rightMostLeafKey.sKey 
				,rightMostLeafKey.keyLength 
				,rightMostLeafKey.pkHash };//��ȻҲ��һ�ַ���ֵ��ʽ
		}
		inline void
			insertAt(uint16_t index, BtreeEntry entry, NodeId childId,
				NodeId rightChild = INVALID_NODEID){
			//Node::insertAtEntryOnly(index, entry);
			//for (uint16_t i = Node::slotuse; i > index;i--) 
			//{
			//	child[i+1] = child[i];
			//}
			//child[index] = childId;
			//if (rightChild != INVALID_NODEID) {
			//	child[index + 1] = rightChild;
			//}
			Node::insertAtEntryOnly(index, entry);
			memmove(&child[index + 1], &child[index],
				sizeof(NodeId)*(slotuse - index));

			child[index] = childId;
			if (rightChild != INVALID_NODEID) {
				child[index + 1] = rightChild;
			}
		}
		/*
		* ���ڵ����һ��entry���ѽ�Ϊ�����ڵ㣬�����м�ڵ�back
		*/
		inline BtreeEntry
			insertSplit(uint16_t insertIndex, BtreeEntry entry,
				InnerNode *emptyRightSibling, NodeId childId,
				NodeId rightChildId = INVALID_NODEID)
		{
			//assert(emptyRightSibling->slotuse == 0);
			//assert(isfull());
			if (emptyRightSibling->slotuse != 0) {
				printf("insertSplit() erorr:Right sibling isn't empty!\n exiting\n");
				exit(0);
			}
			if (!isfull()) {
				printf("insertSplit() erorr:Inner Node isn't full!\n exiting\n");
				exit(0);
			}
			//half=slotuse/2
			uint16_t half = uint16_t(slotuse >> 1);
			if (insertIndex <= half) {
				memmove(&emptyRightSibling->child[0],
					&child[half],
					sizeof(NodeId)*(half + 1));
				//�Ѻ�벿��ת�����Ҳ��ֵܽڵ�
				Node::moveBackEntriesToFrontOf(emptyRightSibling, half);
				insertAt(insertIndex, entry, childId, rightChildId);//����
			}
			else {
				memmove(&emptyRightSibling->child[0],
					&child[half + 1],
					sizeof(NodeId)*(half));
				//��벿���ڲ���ǰת�����Ҳ��ֵܽڵ�
				Node::moveBackEntriesToFrontOf(emptyRightSibling, uint16_t(half - 1));
				emptyRightSibling->insertAt(uint16_t(insertIndex - (half + 1)),
					entry, childId, rightChildId);//����
			}

			// Migrate rightmost leaf key and set a new one for self.
			BtreeEntry back = Node::back();
			if (!rightMostLeafKeyIsInfinite)
				//��ʱrightMostLeafKey��ֵ��δ�ı䣬�����Ҳ��ֵܽڵ㽫�õ���ǰ������Ҷ�Ӽ�
				emptyRightSibling->setRightMostLeafKey(getRightMostLeafKey());

			rightMostLeafKeyIsInfinite = false;
			rightMostLeafKey.sKey = back.sKey;
			rightMostLeafKey.keyLength = back.sKeyLength;
			rightMostLeafKey.pkHash = back.pKeyHash;
			Node::pop_back();
			return back;
		}
		/*
			inbetween�ڵ������ҽڵ���м�ڵ㡣
		*/
		inline void
			mergeIntoRight(InnerNode *right, BtreeEntry inbetween)
		{
			Node::insertAtEntryOnly(slotuse, inbetween);
			//���Ҳ�ڵ��child����ƶ����ڵ�slotuse��λ��
			memmove(&right->child[slotuse], &right->child[0],
				sizeof(NodeId)*(right->slotuse + 1));
			//�����ڵ�child�ƶ����Ҳ�ڵ��λ��
			memmove(&right->child[0], &child[0],
				sizeof(NodeId)*slotuse);
			//��ת�ƽڵ�
			Node::moveBackEntriesToFrontOf(right, slotuse);
			slotuse = 0;
		}
	};
	struct LeafNode : public Node{
		// ˫����ָ��
		NodeId prevleaf, nextleaf;
		LeafNode()
			: Node(0)
			, prevleaf(INVALID_NODEID)
			, nextleaf(INVALID_NODEID)
		{}
		/*
			split���� ���ڵ�ĺ󼸸��ڵ�ת�Ƶ��Ҳ��ֵܽڵ�
		*/
		inline BtreeEntry
			split(uint16_t numEntries, LeafNode *rightSibling) {
			Node::moveBackEntriesToFrontOf(rightSibling, numEntries);
			return Node::back();
		}
		/*
			��ָ���ڵ���뵽ָ��λ��
		*/
		inline void
			insertAt(uint16_t index, BtreeEntry entry) {
			Node::insertAtEntryOnly(index, entry);
		}
		/*
			ɾ��ָ��λ�õĽڵ�
		*/
		inline void
			eraseAt(uint16_t index) {
			Node::eraseAtEntryOnly(index);
		}
		/*
			ƽ�������ڵ�Ľڵ�����half�������ڵ�entry����ƽ��ֵ
			�ұ߱����entry��С��ʱ�����ҽڵ�ת��entry
			��֮��Ȼ
		*/
		inline BtreeEntry
			balanceWithRight(LeafNode *right) {
			uint16_t toMove, half = uint16_t((slotuse + right->slotuse) >> 1);
			//if (right->slotuse > half) {
			//	toMove = uint16_t(right->slotuse - half);
			//	right->moveFrontEntriesToBackOf(this, toMove);
			//}
			//else
			//{
			//	toMove = uint16_t(slotuse- half);
			//	right->moveFrontEntriesToBackOf(this, toMove);
			//}
			if (right->slotuse < slotuse) {
				toMove = uint16_t(half - right->slotuse);
				Node::moveBackEntriesToFrontOf(right, toMove);
			}
			else {
				toMove = uint16_t(half - slotuse);
				right->moveFrontEntriesToBackOf(this, toMove);
			}

			return back();
		}
	};
	std::map<NodeId, Nodeptr> cache;
	Node node1;
	Node node2;
	InnerNode innerNode_1;
	InnerNode innerNode_2;
};

IndexBtree::IndexBtree()
{

}

IndexBtree::~IndexBtree()
{
}
int main() 
{
	BtreeEntry BE_0(1,100,1000);
	BtreeEntry BE_1(1, 101, 1001);
	BtreeEntry BE_2(1, 102, 1002);
	BtreeEntry BE_3(1, 103, 1003);
	IndexBtree btree;
	btree.node1.setAt(0, BE_0);
	btree.node1.setAt(1, BE_1);
	btree.node1.setAt(2, BE_2);
	btree.node1.setAt(3, BE_3);
	getchar();
	btree.node1.moveFrontEntriesToBackOf(&btree.node2, 2);
	//btree.node1.moveBackEntriesToFrontOf(&btree.node2, 2);
	//btree.node.eraseAtEntryOnly(3);
	//btree.node.eraseAtEntryOnly(2);
	btree.node1.eraseAtEntryOnly(1);
	//btree.node.eraseAtEntryOnly(0);
	btree.node1.insertAtEntryOnly(1, BE_1);

	getchar();
}