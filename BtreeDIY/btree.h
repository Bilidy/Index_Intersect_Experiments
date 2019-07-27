#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string> 
#include <map>
#include <iostream>


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
		uint64_t itemcount;// 多少个items
		uint64_t leaves;// 多少叶子节点
		uint64_t innernodes;// 内节点数
		static const uint16_t leafslots = leafslotmax;// 每个叶节点的slot数 leafslotmax = 8
		static const uint16_t innerslots = innerslotmax;// 每个内节点的slot数 innerslotmax = leafslotmax
		// minleafslots = 4;  leafslotmax = leafslots = 8
		// mininnerslots = 4; innerslotmax = innerslots = 8 
		inline tree_stats()
			: itemcount(0),
			leaves(0),
			innernodes(0)
		{}
		inline uint64_t nodes() const {// 返回总node数
			return innernodes + leaves;
		}
		inline double avgfill_leaves() const {// Return the average fill of leaves
			return static_cast<double>((itemcount) / (leaves * leafslots));
		}
	};
	

	//IndexBtree();
	//~IndexBtree();

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
			if (Node::slotuse < n) 
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
			isinnernode() const {
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
			else//index>=Node::slotuse index=4 slotuse=4
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

			uint16_t splitPoint = uint16_t(slotuse - numEntries);

			memmove(&dest->keys[numEntries],
				&dest->keys[0],
				sizeof(KeyInfo)*dest->slotuse);

			memmove(&dest->keys[0],
				&keys[splitPoint],
				sizeof(KeyInfo)*numEntries);

			dest->slotuse = uint16_t(dest->slotuse + numEntries);
			slotuse = uint16_t(slotuse - numEntries);

			/*if (numEntries + dest->slotuse > IndexBtree::innerslotmax) 
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
			}*/
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
			: Node(level)
			, child()
			, rightMostLeafKey()
			, rightMostLeafKeyIsInfinite(true){}
		virtual ~InnerNode() {};

		void 
			setRightMostLeafKey(BtreeEntry entry) {

			rightMostLeafKey.sKey = entry.sKey;
			rightMostLeafKey.keyLength = entry.sKeyLength;
			rightMostLeafKey.pkHash = entry.pKeyHash;

			rightMostLeafKeyIsInfinite = false;
		}
		BtreeEntry 
			getRightMostLeafKey() {
			return { 
				rightMostLeafKey.sKey 
				,rightMostLeafKey.keyLength 
				,rightMostLeafKey.pkHash };//竟然也是一种返回值格式
		}
		inline void
			insertAt(uint16_t index, BtreeEntry entry, NodeId childId,
				NodeId rightChild = INVALID_NODEID){
			Node::insertAtEntryOnly(index, entry);
			memmove(&child[index + 1], &child[index],
				sizeof(NodeId)*(slotuse - index));

			child[index] = childId;
			if (rightChild != INVALID_NODEID) {
				child[index + 1] = rightChild;
			}
		}
		/*
		* 满节点插入一个entry后裂解为两个节点，返回中间节点back
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
				//把后半部分转移至右侧兄弟节点
				Node::moveBackEntriesToFrontOf(emptyRightSibling, half);
				insertAt(insertIndex, entry, childId, rightChildId);//插入
			}
			else {
				memmove(&emptyRightSibling->child[0],
					&child[half + 1],
					sizeof(NodeId)*(half));
				//后半部分在插入前转移至右侧兄弟节点
				Node::moveBackEntriesToFrontOf(emptyRightSibling, uint16_t(half - 1));
				emptyRightSibling->insertAt(uint16_t(insertIndex - (half + 1)),
					entry, childId, rightChildId);//插入
			}

			// Migrate rightmost leaf key and set a new one for self.
			BtreeEntry back = Node::back();
			if (!rightMostLeafKeyIsInfinite)
				//此时rightMostLeafKey的值尚未改变，所以右侧兄弟节点将得到先前的最右叶子键
				emptyRightSibling->setRightMostLeafKey(getRightMostLeafKey());

			rightMostLeafKeyIsInfinite = false;
			rightMostLeafKey.sKey = back.sKey;
			rightMostLeafKey.keyLength = back.sKeyLength;
			rightMostLeafKey.pkHash = back.pKeyHash;
			Node::pop_back();
			return back;
		}
		inline void
			eraseAt(uint16_t index)
		{
			// Special case for inner nodes. Since inner nodes can have
			// n + 1 keys and n + 1 pointers, erasing the last pointer is
			// equivalent to promoting the nth pointer to the last pointer.
			if (index == slotuse) {
				slotuse--;
				rightMostLeafKey.sKey = keys[slotuse].sKey;
				rightMostLeafKey.keyLength = keys[slotuse].keyLength;
				rightMostLeafKey.pkHash = keys[slotuse].pkHash;
			}
			else {
				Node::eraseAtEntryOnly(index);
				memmove(&child[index], &child[index + 1],
					sizeof(NodeId)*(slotuse + 1 - index));
			}
		}
		inline NodeId
			getChildAt(uint16_t index) const{
			if (index > Node::slotuse) {
				printf("Attempted to get a child at an index larger than slotuse");
			}

			return child[index];
		}
		inline BtreeEntry
			balanceWithRight(InnerNode *right, BtreeEntry inbetween)
		{
			uint16_t entriesToMove, childrenToMove;
			uint16_t nLeft = slotuse;
			uint16_t nRight = right->slotuse;
			uint16_t half = uint16_t((nLeft + nRight) >> 1);

			// insert the missing link
			uint16_t extraRightChild;
			if (isfull()) {
				right->Node::insertAtEntryOnly(0, inbetween);
				extraRightChild = 1;
			}
			else {
				Node::insertAtEntryOnly(slotuse, inbetween);
				extraRightChild = 0;
			}

			// Move children and perform balance
			if (nRight < nLeft) {
				entriesToMove = uint16_t(half - nRight);
				childrenToMove = uint16_t(entriesToMove + extraRightChild);

				memmove(&right->child[childrenToMove], &right->child[0],
					sizeof(NodeId)*(nRight + 1));
				memmove(&right->child[0], &child[nLeft - childrenToMove + 1],
					sizeof(NodeId)*(childrenToMove));

				Node::moveBackEntriesToFrontOf(right, entriesToMove);
			}
			else {
				entriesToMove = uint16_t(half - nLeft);

				memmove(&child[nLeft + 1], &right->child[0],
					sizeof(NodeId)*(entriesToMove));

				memmove(&right->child[0], &right->child[entriesToMove],
					sizeof(NodeId)*(nRight - entriesToMove + 1));

				right->moveFrontEntriesToBackOf(this, entriesToMove);
			}

			BtreeEntry back = Node::back();
			rightMostLeafKey.keyLength = back.sKeyLength;
			rightMostLeafKey.pkHash = back.pKeyHash;
			Node::pop_back();
			return back;
		}
		/*
			inbetween节点是左右节点的中间节点。
		*/
		inline void
			mergeIntoRight(InnerNode *right, BtreeEntry inbetween)
		{
			Node::insertAtEntryOnly(slotuse, inbetween);
			//把右侧节点的child向后移动左侧节点slotuse个位置
			memmove(&right->child[slotuse], &right->child[0],
				sizeof(NodeId)*(right->slotuse + 1));
			//将左侧节点child移动到右侧节点空位中
			memmove(&right->child[0], &child[0],
				sizeof(NodeId)*slotuse);
			//将转移节点
			Node::moveBackEntriesToFrontOf(right, slotuse);
			slotuse = 0;
		}
		
};
	struct LeafNode : public Node{
		// 双链表指针
		NodeId prevleaf, nextleaf;
		LeafNode()
			: Node(0)
			, prevleaf(INVALID_NODEID)
			, nextleaf(INVALID_NODEID)
		{}
		virtual ~LeafNode() {};
		//inline LeafNode& operator=(const LeafNode&other)
		//{
		//	level = other.level;
		//	slotuse = other.slotuse;
		//	keyStorageUsed = other.keyStorageUsed;
		//	prevleaf = other.prevleaf;
		//	nextleaf = other.nextleaf;
		//	return *this;
		//}
		/*
			split函数 将节点的后几个节点转移到右侧兄弟节点
		*/
		inline BtreeEntry split(uint16_t numEntries, LeafNode *rightSibling) {
			Node::moveBackEntriesToFrontOf(rightSibling, numEntries);
			return Node::back();
		}
		/*
			将指定节点插入到指定位置
		*/
		inline void
			insertAt(uint16_t index, BtreeEntry entry) {
			Node::insertAtEntryOnly(index, entry);
		}
		/*
			删除指定位置的节点
		*/
		inline void
			eraseAt(uint16_t index) {
			Node::eraseAtEntryOnly(index);
		}
		/*
			平衡两个节点的节点数，half是两个节点entry数的平均值
			右边比左边entry数小的时候，向右节点转移entry
			反之亦然
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
		inline void
		mergeIntoRight(LeafNode *right) {
			Node::moveBackEntriesToFrontOf(right,slotuse);
			slotuse = 0;
		}

	};
	std::map<NodeId, Nodeptr> nodeCache;
	explicit inline IndexBtree(uint64_t tableId)
		: m_stats()
		, treeTableId(tableId)
		, nextNodeId(ROOT_ID)
		, m_rootId(ROOT_ID)
		, numEntries(0)
		, nodeCache()
	{ }
	explicit inline IndexBtree(uint64_t tableId,uint64_t nextNodeId)
		: m_stats()
		, treeTableId(tableId)
		, nextNodeId(nextNodeId)
		, m_rootId(ROOT_ID)
		, numEntries(0)
		, nodeCache()
	{ }
	inline ~IndexBtree() {}
	/*static bool
		isGreaterOrEqual(Buffer* nodeObjectValue, BtreeEntry compareEntry) {}*/
public:
	NodeId
		getNextNodeId() {
		return nextNodeId;
	}
	void
		setNextNodeId(NodeId newNodeId) {
		nextNodeId = newNodeId;
	}
public:
	void
	clear_fast() {
		if (nextNodeId > ROOT_ID) {
			
			nextNodeId = ROOT_ID;
			m_stats = tree_stats();
			
			nodeCache.clear();
		}
	}
	/**
	* Destroys the B+ tree by recursively freeing all the nodes *expensive*
	*/
	//void
	//	clear() {
	//	if (nextNodeId > ROOT_ID) {
	//		clear_recursive(ROOT_ID);
	//		m_stats = tree_stats();
	//	}

	//	flush();
	//}
	//void
	//clear_recursive(NodeId nodeId) {
	//	Node *n = readNode(nodeId);

	//	if (n->isinnernode()) {
	//		const InnerNode *innernode = static_cast<const InnerNode*>(n);
	//		for (uint16_t slot = 0; slot <= innernode->slotuse; ++slot)
	//			clear_recursive(innernode->getChildAt(slot));
	//	}

	//	freeNode(nodeId);
	//}
	void clear(){
		if (nextNodeId > ROOT_ID) {
			clear_recursive(ROOT_ID);
			m_stats = tree_stats();
		}
	}

	class iterator;//前置声明

	/*返回指向B+Tree最左侧值的迭代器*/
	iterator begin() {
		if (nextNodeId <= ROOT_ID)
			return iterator(this, INVALID_NODEID, 0);
		NodeId currentId = m_rootId;
		Node *n = readNode(currentId);
		while (!n->isLeaf()) {
			const InnerNode *inner = static_cast<const InnerNode*>(n);

			currentId = inner->getChildAt(0);
			n = readNode(currentId);
		}
		return iterator(this, currentId, 0);
	}
	/*返回指向B+Tree最右侧节点第一个不可用槽的迭代器*/
	iterator end() {
		return iterator(this, INVALID_NODEID, 0);
	}
	uint64_t size() const {
		if (nextNodeId <= ROOT_ID)
			return 0;
		else
			return m_stats.itemcount;
	}
	bool empty() const {
		return (nextNodeId == ROOT_ID);
	}
	bool exists(const BtreeEntry &key) {
		return find(key) != end();
	}
	/*返回指向B+树指定值位置的迭代器*/
	iterator find(const BtreeEntry &key) {
		if (nextNodeId <= ROOT_ID)
			return end();
		NodeId currId = m_rootId;
		Node *n = readNode(m_rootId);
		while (!n->isLeaf())
		{
			const InnerNode *inner= static_cast<const InnerNode*>(n);
			uint16_t slot= findEntryGE(inner, key);

			currId = inner->getChildAt(slot);
			n = readNode(currId);
		}
		uint16_t slot = findEntryGE(n, key);
		return (slot<n->slotuse&&key_equal(n->getAt(slot),key))
			? iterator(this, currId, slot) :end();
	}
	iterator lower_bound(const BtreeEntry key) {
		if (nextNodeId <= ROOT_ID)
			return end();

		Node *n = readNode(m_rootId);
		NodeId childId = m_rootId;
		while (!n->isLeaf()) {
			const InnerNode *inner = static_cast<const InnerNode*>(n);
			uint16_t slot = findEntryGE(inner, key);
			childId = inner->getChildAt(slot);
			n = readNode(childId);
		}

		const LeafNode *leaf= static_cast<const LeafNode*>(n);
		uint16_t slot = findEntryGE(leaf, key);

		if (slot >= leaf->slotuse)
			return end();
		else
			return iterator(this, childId, slot);
	}
	iterator upper_bound(const BtreeEntry& key) {
		if (nextNodeId <= ROOT_ID)
			return end();

		Node *n = readNode(m_rootId);
		NodeId childId = m_rootId;
		while (!n->isLeaf()) {
			InnerNode *inner = static_cast<InnerNode*>(n);
			uint16_t slot = findEntryGreater(inner, key);
			childId = inner->getChildAt(slot);
			n = readNode(childId);
		}

		const LeafNode *leaf = static_cast<const LeafNode*>(n);
		uint16_t slot = findEntryGreater(leaf, key);

		if (slot >= leaf->slotuse)
			return end();
		else
			return iterator(this, childId, slot);
	}
	/*key键的entr存在数量*/
	uint64_t count(const BtreeEntry &key) const {
		if (nextNodeId <= ROOT_ID)
			return 0;

		NodeId childId = m_rootId;
		Node *n = readNode(m_rootId); 
		
		while (!n->isLeaf()){
			const InnerNode *inner = static_cast<InnerNode *>(n);
			uint16_t slot = findEntryGE(inner, key);
			childId = inner->getChildAt(slot);
			n = readNode(childId);
		}

		const LeafNode *leaf = static_cast<const LeafNode*>(n);
		uint16_t slot = findEntryGE(leaf,key);
		uint64_t counter = 0;

		NodeId currentLeafId = childId;
		while (currentLeafId != INVALID_NODEID 
			&& slot< leaf->slotuse
			&& key_equal(key,leaf->getAt(slot))) {
			counter++;
			if (++slot=leaf->slotuse)
			{
				slot = 0;
				currentLeafId = leaf->nextleaf;
				if (currentLeafId == INVALID_NODEID) {
					break;
				}
				leaf= static_cast<const LeafNode*>(readNode(currentLeafId));
			}
		}
		return counter;
	}
	void insert(const BtreeEntry entry) {
		if (nextNodeId == ROOT_ID) {
			LeafNode rootemp;
			LeafNode *root = (LeafNode*)malloc(sizeof(LeafNode)); //rootBuffer.emplaceAppend<LeafNode>(&rootBuffer);
			*root = rootemp;
			root->insertAt(0, entry);
			writeNode(root, ROOT_ID);
			nextNodeId = ROOT_ID + 1;
			m_stats.leaves = 1;
		}
		else {
			ChildUpdateInfo info;
			insertDescend(ROOT_ID, entry, &info);

			// Root node was split
			if (info.childSplit) {
				uint16_t rootLevel = uint16_t(info.getChildLevel() + 1);

				InnerNode newRoottemp(rootLevel);
				InnerNode *newRoot = (InnerNode*)malloc(sizeof(InnerNode)); //rootBuffer.emplaceAppend<InnerNode>(&rootBuffer, rootLevel);
				*newRoot = newRoottemp;

				newRoot->insertAt(0,
					info.newChild,
					info.newChildId,
					info.rightSiblingId);
				writeNode(newRoot, ROOT_ID);
				m_stats.innernodes++;
			}
		}
		m_stats.itemcount++;
	}
	bool erase(BtreeEntry entry) {
		//if (selfverify) verify();

		// The tree is empty; do nothing.
		if (nextNodeId <= ROOT_ID)
			return false;

		EraseUpdateInfo info;
		bool success = eraseOneDescend(entry, m_rootId, NULL, 0, &info);

		//flush();
		//if (selfverify) verify();
		return success;
	}
private:
	inline bool key_less(const BtreeEntry a, const BtreeEntry b) const {
		//a<b
		return ((a.sKey - b.sKey) == 0 )
			?(a.pKeyHash<b.pKeyHash):(a.sKey < b.sKey);
	}
	static bool key_less_static(const BtreeEntry a, const BtreeEntry b) {
		return ((a.sKey - b.sKey) == 0)
			? (a.pKeyHash < b.pKeyHash) : ((a.sKey < b.sKey));
	}
	inline bool key_lessequal(const BtreeEntry a, const BtreeEntry b) const {
		return !key_less(b, a);
	}
	inline bool key_greater(const BtreeEntry a, const BtreeEntry &b) const {
		return key_less(b, a);
	}
	inline bool key_greaterequal(const BtreeEntry a, const BtreeEntry b) const {
		return !key_less(a, b);
	}
	static bool key_greaterequal_static(const BtreeEntry a, const BtreeEntry b) {
		return !key_less_static(a, b);
	}
	inline bool key_equal(const BtreeEntry a, const BtreeEntry b) const {
		return a.sKey==b.sKey&& a.pKeyHash == b.pKeyHash;
	}
	inline uint16_t findEntryGE(const Node *n, BtreeEntry entry) const{
		if (useBinarySearch) {
			if (n->slotuse == 0) {
				return 0;
			}
			uint16_t lo = 0, hi = n->slotuse;
			while (lo < hi){
				uint16_t mid = uint16_t((lo + hi) >> 1);
				if (key_lessequal(entry,n->getAt(mid))){
					hi = mid;
				}
				else {
					lo = uint16_t(mid+1);
				}
			}
			if (selfverify)
			{
				uint16_t i = 0;
				while (i < n->slotuse&&key_less(n->getAt(i), entry))
				{
					i++;
				}
				if (i != lo) {
					std::cout<<"findEntryGE selfverify not equal"<<std::endl;
				}
			}
			return lo;
		} else {
			uint16_t lo = 0;
			while (lo < n->slotuse && key_less(n->getAt(lo), entry)) ++lo;
			return lo;
		}
	}
	inline uint16_t findEntryGreater(const Node *n, const BtreeEntry entry) const {
		if (useBinarySearch) {
			if (n->slotuse == 0)
				return 0;

			uint16_t lo = 0, hi = n->slotuse;
			while (lo < hi) {
				uint16_t mid = uint16_t((lo + hi) >> 1);
				if (key_less(entry, n->getAt(mid))) {
					hi = mid; // key < mid
				} else {
					lo = uint16_t(mid + 1); // key >= mid
				}
			}
			if (selfverify)
			{
				uint16_t i = 0;
				while (i < n->slotuse && key_lessequal(n->getAt(i), entry)) ++i;
				if (i != hi) {
					std::cout << "findEntryGreater selfverify not equal" << std::endl;
				}
			}
			return lo;
		}
		else {
			uint16_t lo = 0;
			while (lo < n->slotuse && key_lessequal(n->getAt(lo), entry)) ++lo;
			return lo;
		}
	}

	inline Node*
		readNode(NodeId nodeId) const {
		Nodeptr n=NULL;
		std::map<NodeId, Nodeptr>::const_iterator iter;
		iter=nodeCache.find(nodeId);
		if (iter != nodeCache.end()) {
			//std::cout << "Find, the value is " << iter->second << std::endl;
			n = iter->second;
		}
		else
			std::cout << "readNode error:Do not Find NodeId:"<< nodeId << std::endl;
		return n;
	}
	inline NodeId
		writeNode(Node *node, NodeId nodeId = INVALID_NODEID) {
		if (nodeId == INVALID_NODEID)
			nodeId = nextNodeId++;
		Node* nodeptr = node;
		nodeCache[nodeId] = node;
		return nodeId;
	}
	void
		clear_recursive(NodeId nodeId) {
		Node *n = readNode(nodeId);
		if (n->isinnernode()) {
			const InnerNode *innernode = static_cast<const InnerNode*>(n);
			for (uint16_t slot = 0; slot <= innernode->slotuse; ++slot)
				clear_recursive(innernode->getChildAt(slot));
		}
		freeNode(nodeId);
	}
	inline void
		freeNode(NodeId nodeId){
		Nodeptr ptr = nodeCache[nodeId];
		if(ptr)
			free(ptr);

		if (nodeCache.erase(nodeId)){
		}
		if (nodeId == m_rootId)
			nextNodeId = ROOT_ID;
	}
private:
	struct ChildUpdateInfo {
		BtreeEntry newChild;	//父节点用以索引孩子的临时节点
		NodeId newChildId;		//孩子节点的NodeId
		NodeId rightSiblingId;	//孩子节点的右兄弟
		uint16_t childLevel;	//孩子结点的层次
		bool childSplit;		//孩子是否发生分裂
		bool rightMostKeyUpdated;//表明最右侧叶节点被更底层更新了，需要向上传播给父节点
		BtreeEntry rightMostLeafKey;
		ChildUpdateInfo()
			: newChild()
			, newChildId(INVALID_NODEID)
			, rightSiblingId(INVALID_NODEID)
			, childLevel(0)
			, childSplit(false)
			, rightMostKeyUpdated(false)
			, rightMostLeafKey() {}
		void setSplit(BtreeEntry currLastEntry, NodeId currChildId,
			NodeId rightSiblingId, uint16_t level) {

			newChild.sKeyLength = currLastEntry.sKeyLength;
			newChild.pKeyHash = currLastEntry.pKeyHash;
			newChild.sKey = currLastEntry.sKey;

			newChildId = currChildId;
			this->rightSiblingId = rightSiblingId;
			childLevel = level;
			childSplit = true;

		}
		void setLastKeyUpdated(BtreeEntry lastKeyIn) {

			rightMostKeyUpdated = true;
			rightMostLeafKey.sKey = lastKeyIn.sKey;
			rightMostLeafKey.sKeyLength = lastKeyIn.sKeyLength;
			rightMostLeafKey.pKeyHash = lastKeyIn.pKeyHash;

		}
		void clear() {
			childSplit = false;
		}
		BtreeEntry getChild() {
			return newChild;
		}
		NodeId getRightSiblingId() {
			return rightSiblingId;
		}
		uint16_t getChildLevel() {
			return childLevel;
		}
	};
	struct EraseUpdateInfo {
		enum updateOperation {
			setLeft,
			setCurr,
			delLeft,
			delCurr,
			noOp
		};
		bool lastEntryUpdated;
		BtreeEntry lastEntry;
		enum updateOperation op;
		BtreeEntry setEntry;
		EraseUpdateInfo()
			: lastEntryUpdated(false)
			, lastEntry()
			, op(noOp)
			, setEntry()
		{}
		void setSetEntry(BtreeEntry e) {
			copyKeyInternal(&setEntry, &e);
		}
		void setLastKey(BtreeEntry e) {
			lastEntryUpdated = true;
			copyKeyInternal(&lastEntry, &e);
		}
		void copyKeyInternal(BtreeEntry *to, BtreeEntry *from) {
			*to = *from;
		}
		void clearOp() {
			op = noOp;
		}
	};
	void insertDescend(NodeId currentId, BtreeEntry entry, ChildUpdateInfo *updateInfo) 
	{
		Node *n = readNode(currentId);

		uint16_t insertIndex=findEntryGE(n,entry);
		if (n->isinnernode()) {
			InnerNode* inner = (InnerNode*)n;
			updateInfo->clear();
			insertDescend(inner->getChildAt(insertIndex), entry, updateInfo);
			bool rightMostKeyUpdated = false;
			if (updateInfo->rightMostKeyUpdated && !inner->rightMostLeafKeyIsInfinite) {
				
				if (insertIndex < inner->slotuse)
				{
					//插入位置并不是在node的最后slot，并没有改变上层索引大小，停止向上传播；
					updateInfo->rightMostKeyUpdated = false;
					//?????
					if (key_less(inner->getAt(insertIndex), updateInfo->rightMostLeafKey))
					{
						rightMostKeyUpdated = true;
						inner->setAt(insertIndex, updateInfo->rightMostLeafKey);
					}
				}
				else 
				if (key_less(inner->getRightMostLeafKey(),
					updateInfo->rightMostLeafKey)) {
					rightMostKeyUpdated = true;
					inner->setRightMostLeafKey(updateInfo->rightMostLeafKey);
				}
			}
			//子节点没有分裂操作，很棒，结束了。
			if (!updateInfo->childSplit) {
				if (rightMostKeyUpdated)
					writeNode(inner, currentId);

				return;
			}
			//哎~!有分裂操作，而且inner没满，要把下层分裂的插入到inner中
			else 
			{
				if (!inner->isfull())
				{
					inner->insertAt(insertIndex,
						updateInfo->newChild,
						updateInfo->newChildId,
						updateInfo->rightSiblingId);
					writeNode(inner, currentId);
					updateInfo->clear();
					return;
				}
				// Inner is full, need to split;
				// 首先申请新的空间作为分裂后右侧兄弟
				//?????这里可能有问题；
				//InnerNode newRightSiblingtmp;
				InnerNode newRightSiblingtemp(inner->level);
				InnerNode *newRightSibling = (InnerNode*)malloc(sizeof(InnerNode));
				*newRightSibling = newRightSiblingtemp;
				newRightSibling->level = inner->level;
				m_stats.innernodes++;
				//做满节点插入（要分裂的）
				BtreeEntry newLastEntry = inner->insertSplit(insertIndex,
					updateInfo->getChild(),
					newRightSibling,
					updateInfo->newChildId,
					updateInfo->getRightSiblingId());
				//如果根节点也被分裂了那就使用一个新的NodeId重新当前根
				NodeId rightId = (currentId == m_rootId) ? nextNodeId++ : currentId;
				NodeId leftId = nextNodeId++;

				updateInfo->setSplit(newLastEntry, leftId, rightId, inner->level);
				if (!newRightSibling->rightMostLeafKeyIsInfinite)
					updateInfo->setLastKeyUpdated(newRightSibling->getRightMostLeafKey());
				//存入节点
				writeNode(inner, leftId);
				writeNode(newRightSibling, rightId);
			}
		}
		else//n不是内部节点，是叶子诶；
		{
			LeafNode *leaf = static_cast<LeafNode*>(n);

			//插入位置在节点最右
			if(insertIndex==leaf->slotuse)
				updateInfo->setLastKeyUpdated(entry);
			//不需要分裂，直接插入，插入后退出
			if (!leaf->isfull()) 
			{
				leaf->insertAt(insertIndex, entry);
				writeNode(leaf, currentId);
				updateInfo->clear();

				return;
			}
			else //叶子是满的，需要分裂。从中间分裂；
			{
				uint16_t midpoint = uint16_t(leaf->slotuse >> 1);
				LeafNode newRightSiblingtemp;
				LeafNode *newRightSibling = (LeafNode*)malloc(sizeof(LeafNode));
				*newRightSibling = newRightSiblingtemp;

				m_stats.leaves++;
				leaf->split(midpoint, newRightSibling);
				if (insertIndex<leaf->slotuse){
					leaf->insertAt(insertIndex, entry);
				}
				else{
					newRightSibling->insertAt(uint16_t(insertIndex - midpoint), entry);
				}

				NodeId rightId = (currentId == m_rootId) ? nextNodeId++ : currentId;
				NodeId leftId = nextNodeId++;

				updateInfo->setSplit(leaf->back(), leftId, rightId, leaf->level);
				updateInfo->setLastKeyUpdated(newRightSibling->back());
				//调整叶子节点的nextleaf和prevleaf
				newRightSibling->nextleaf = leaf->nextleaf;
				newRightSibling->prevleaf = leftId;
				leaf->nextleaf = rightId;

				writeNode(leaf, leftId);
				writeNode(newRightSibling, rightId);

				//更新左兄弟
				if (leaf->prevleaf != INVALID_NODEID) {
					LeafNode *prevLeaf = static_cast<LeafNode*>(readNode(leaf->prevleaf));
					prevLeaf->nextleaf = leftId;
					writeNode(prevLeaf, leaf->prevleaf);
				}
			}
			
		}
	}
	bool eraseOneDescend(BtreeEntry key, NodeId currentId,
			const InnerNode *parent, uint16_t parentSlot,
			EraseUpdateInfo *info) {
		Node *node = readNode(currentId);
		uint16_t slot = findEntryGE(node, key);//在孩子结点中所属槽位
		bool dirty = false;

		if (node->isLeaf())
		{
			if (slot >= node->slotuse || !key_equal(key, node->getAt(slot)))
				return false;
			LeafNode *currLeaf = static_cast<LeafNode*>(node);
			currLeaf->eraseAt(slot);
			m_stats.itemcount--;
			dirty = true;

			if (slot == currLeaf->slotuse && currentId != ROOT_ID)
				info->setLastKey(currLeaf->back());
		}
		else
		{
			InnerNode *inner = static_cast<InnerNode*>(node);
			NodeId childId = inner->getChildAt(slot);//得到孩子所在槽位指向的孩子节点id

			info->clearOp();
			bool success = eraseOneDescend(key, childId, inner, slot,info);
			if (!success){
				return false;
			}
			if (info->lastEntryUpdated && !inner->rightMostLeafKeyIsInfinite) {
				if (slot < inner->slotuse) {
					info->lastEntryUpdated = false;
					if (key_less(inner->getAt(slot),info->lastEntry))
					{
						inner->setAt(slot, info->lastEntry);
						dirty = true;
					}
				}
				else if (key_less(inner->getRightMostLeafKey(),info->lastEntry)) 
				{
					inner->setRightMostLeafKey(info->lastEntry);
					dirty = true;
				}
			}
			dirty |= (info->op < EraseUpdateInfo::noOp);
			uint16_t leftSlot = uint16_t(slot - 1);
			switch (info->op) {
				case EraseUpdateInfo::setLeft:
					inner->setAt(leftSlot, info->setEntry);
					break;
				case EraseUpdateInfo::setCurr:
					inner->setAt(slot, info->setEntry);
					break;
				case EraseUpdateInfo::delLeft:
					inner->eraseAt(leftSlot);
					break;
				case EraseUpdateInfo::delCurr:
					inner->eraseAt(slot);
					break;
				default:
					break;
			}
			info->clearOp();
		}
		
		if (node->isunderflow() && !(currentId == m_rootId && node->slotuse >= 1))
		{
			handleUnderflowAndWrite(node, currentId, parent, parentSlot, info);
		}
		else if (true){
			writeNode(node, currentId);
		}
		return true;
	}
	inline void
		handleUnderflowAndWrite(Node *curr, NodeId currId,
			const InnerNode *parent, uint16_t parentSlot,
			EraseUpdateInfo *info){
		//存在一种特殊情况，作为currId为根节点时，我们允许存在下溢
		if (currId == m_rootId && curr->slotuse == 0) {
			if (curr->isLeaf()) {
				freeNode(m_rootId);
				m_stats.leaves--;
				return;
			}
			//根节点并不是叶子节点，所以我们需要提升最后一个孩子节点
			InnerNode *inner = static_cast<InnerNode*>(curr);
			NodeId childId = inner->getChildAt(0);
			std::map<NodeId, Nodeptr>::iterator it = nodeCache.find(childId);
			Node *newRoot;
			if (it == nodeCache.end()) {
				newRoot = readNode(childId);
			}
			else {
				newRoot = static_cast<Node*>(
					it->second);

				// Readjust buffer????
				//newRoot->reinitFromRead(&logBuffer, it->second);
			}
			writeNode(newRoot, m_rootId);
			freeNode(childId);
			m_stats.innernodes--;
			return;
		}
		uint16_t rightSlot = uint16_t(parentSlot + 1);
		uint16_t leftSlot = uint16_t(parentSlot - 1);

		/// Case 1: There is only the right sibling to merge/balance with
		// Case 1:只存在右侧兄弟合并/平衡
		if (parentSlot == 0) {
			NodeId rightId = parent->getChildAt(rightSlot);
			Node *right = readNode(rightId);

			// Edge case where current key was updated in the current delete op
			BtreeEntry currParentKey = (info->lastEntryUpdated)
				? info->lastEntry : parent->getAt(parentSlot);

			if (right->isfew()) {
				merge(curr, right, currParentKey);
				info->op = info->delCurr;
				writeNode(right, rightId);
				freeNode(currId);
			}
			else {
				BtreeEntry newLastKey = balance(curr, right, currParentKey);
				info->op = info->setCurr;
				info->setSetEntry(newLastKey);
				writeNode(curr, currId);
				writeNode(right, rightId);
			}
		}
		// Case 2: There is only left sibling to merge/balance with.
		// Case 2:只有左侧兄弟结点合并/平衡
		else if (parentSlot == parent->slotuse) {
			NodeId leftId = parent->getChildAt(leftSlot);
			Node *left = readNode(leftId);
			BtreeEntry leftParentKey = parent->getAt(leftSlot);

			if (left->isfew()) {
				merge(left, curr, leftParentKey);
				info->op = info->delLeft;
				writeNode(curr, currId);
				freeNode(leftId);
			}
			else {
				BtreeEntry newLastKey = balance(left, curr, leftParentKey);
				info->op = info->setLeft;
				info->setSetEntry(newLastKey);
				writeNode(left, leftId);
				writeNode(curr, currId);
			}
		}
		// Case 3: Both siblings exist
		// Case 3: 两个兄弟都存在
		else {
			NodeId leftId = parent->getChildAt(leftSlot);
			NodeId rightId = parent->getChildAt(rightSlot);
			Node *left = readNode(leftId);
			Node *right = readNode(rightId);
			BtreeEntry leftParentKey = parent->getAt(leftSlot);

			// Edge case where current key was updated in the current delete op
			//
			BtreeEntry currParentKey =
				(info->lastEntryUpdated) ? info->lastEntry : parent->getAt(parentSlot);

			// Case 3a: Right is about to underflow so merge with it.
			// Case 3a: 右侧结点将要发生下溢出。
			if (right->isfew()) {
				merge(curr, right, currParentKey);
				info->op = info->delCurr;
				writeNode(right, rightId);
				freeNode(currId);
			}
			// Case 3b: Left is about to underflow so merge with it
			else if (left->isfew()) {
				merge(left, curr, leftParentKey);
				info->op = info->delLeft;
				writeNode(curr, currId);
				freeNode(leftId);
			}
			// Case 3c: Both left and right has extras, so balance with the
			// one with more extras.
			else if (left->slotuse >= right->slotuse) {
				BtreeEntry newLastKey = balance(left, curr, leftParentKey);
				info->op = info->setLeft;
				info->setSetEntry(newLastKey);
				writeNode(left, leftId);
				writeNode(curr, currId);
			}
			else {
				BtreeEntry newLastKey = balance(curr, right, currParentKey);
				info->op = info->setCurr;
				info->setSetEntry(newLastKey);
				writeNode(right, rightId);
				writeNode(curr, currId);
			}
		}

		// If merge occurred, decrement the node count
		if (info->op == info->delCurr || info->op == info->delLeft) {
			if (curr->level == 0) {
				m_stats.leaves--;
			}
			else {
				m_stats.innernodes--;
			}
		}
	}
	inline BtreeEntry
		balance(Node *left, Node *right, BtreeEntry leftParentEntry)
	{
		if (left->isinnernode()) {
			InnerNode *leftIn = static_cast<InnerNode*>(left);
			InnerNode *rightIn = static_cast<InnerNode*>(right);
			return leftIn->balanceWithRight(rightIn, leftParentEntry);
		}
		else {
			LeafNode *leafLeaf = static_cast<LeafNode*>(left);
			LeafNode *rightLeaf = static_cast<LeafNode*>(right);
			return leafLeaf->balanceWithRight(rightLeaf);
		}
	}
	inline void
		merge(Node *left, Node *right, BtreeEntry leftParentEntry)
	{

		// If these are leaf nodes, we need to adjust pointers
		if (left->isLeaf()) {
			LeafNode *leftLeaf = static_cast<LeafNode*>(left);
			LeafNode *rightLeaf = static_cast<LeafNode*>(right);
			leftLeaf->mergeIntoRight(rightLeaf);

			rightLeaf->prevleaf = leftLeaf->prevleaf;
			if (leftLeaf->prevleaf != INVALID_NODEID) {
				LeafNode *prev = static_cast<LeafNode*>(
					readNode(leftLeaf->prevleaf));
				prev->nextleaf = leftLeaf->nextleaf;
				writeNode(prev, leftLeaf->prevleaf);
			}
		}
		else {
			InnerNode *leftInner = static_cast<InnerNode*>(left);
			InnerNode *rightInner = static_cast<InnerNode*>(right);
			leftInner->mergeIntoRight(rightInner, leftParentEntry);
		}
	}
public:
	class iterator {
		typedef std::bidirectional_iterator_tag iterator_category;

	private:
		
		IndexBtree* parentBtree;	//当前指向的B树指针
		NodeId currentNodeId;		//当前引用的叶子节点NodeId
		uint16_t currslot;			//当前的slot
		const LeafNode* currnode;	//当前叶子节点的指针
		BtreeEntry tempEntry;		//临时BtreeEntry

	public:

		inline iterator(IndexBtree* tree = NULL)
			: parentBtree(tree)
			, currentNodeId(INVALID_NODEID)
			, currslot(0)
			, currnode(NULL)
			, tempEntry()
		{ }
		inline iterator(IndexBtree *tree, NodeId nodeId, uint16_t slot = 0)
			: parentBtree(tree)
			, currentNodeId(nodeId)
			, currslot(slot)
			, currnode(NULL)
			, tempEntry()
		{ }
		inline iterator(const iterator &it)
			: parentBtree(it.parentBtree)
			, currentNodeId(it.currentNodeId)
			, currslot(it.currslot)
			, currnode(NULL)
			, tempEntry()
		{ }
		inline iterator&
			operator=(const iterator &it){
			parentBtree = it.parentBtree;
			currentNodeId = it.currentNodeId;
			currslot = it.currslot;
			currnode = it.currnode;
			return *this;
		}
		inline ~iterator() { }
		/**/
		inline BtreeEntry&
			operator*(){
			if (!currnode){
				currnode= static_cast<const LeafNode *>(parentBtree->readNode(currentNodeId));
			}
			tempEntry = currnode->getAt(currslot);
		}
		/**/
		inline BtreeEntry*
			operator->(){
			if (!currnode) {
				currnode = static_cast<const LeafNode *>(parentBtree->readNode(currentNodeId));
			}
			tempEntry = currnode->getAt(currslot);
			return &tempEntry;
		}
		inline iterator&
			operator++(){
			if (!currnode)
				currnode = static_cast<const LeafNode*>(
					parentBtree->readNode(currentNodeId));

			if (currslot + 1< currnode->slotuse){
				currslot++;
			}
			else if (currnode->nextleaf!=INVALID_NODEID) {
				
				currentNodeId = currnode->nextleaf;
				currnode = NULL;
				currslot = 0;
			}
			else{
				currentNodeId = INVALID_NODEID;
				currslot = 0;
				currnode = NULL;
			}
			return *this;
		}
		inline iterator
			operator++(int){
			iterator tmp = *this;   // copy ourselves

			if (!currnode)
				currnode = static_cast<const LeafNode*>(
					parentBtree->readNode(currentNodeId));

			if (currslot + 1 < currnode->slotuse) {
				currslot++;
			}
			else if (currnode->nextleaf != INVALID_NODEID) {

				currentNodeId = currnode->nextleaf;
				currnode = NULL;
				currslot = 0;
			}
			else {
				currentNodeId = INVALID_NODEID;
				currslot = 0;
				currnode = NULL;
			}
			return tmp;
		}
		inline iterator&
			operator--(){
			if (!currnode){
				currnode = static_cast<const LeafNode*>(
					parentBtree->readNode(currentNodeId));
			}
			if (currslot > 0){
				currslot--;
			}
			else if (currnode->prevleaf!=INVALID_NODEID) {
				currentNodeId = currnode->prevleaf;
				currslot= uint16_t(currnode->slotuse - 1);
				currnode = NULL;
			}
			else{
				currslot = 0;
				currnode = NULL;
			}
			return *this;
		}
		inline iterator
			operator--(int)
		{
			iterator tmp = *this;

			if (!currnode) {
				currnode = static_cast<const LeafNode*>(
					parentBtree->readNode(currentNodeId));
			}
			if (currslot > 0) {
				currslot--;
			}
			else if (currnode->prevleaf != INVALID_NODEID) {
				currentNodeId = currnode->prevleaf;
				currslot = uint16_t(currnode->slotuse - 1);
				currnode = NULL;
			}
			else {
				currslot = 0;
				currnode = NULL;
			}
			return tmp;
		}
		inline bool
			operator==(const iterator& x) const{
			return (x.parentBtree == parentBtree)
				&& (x.currentNodeId == currentNodeId)
				&& (x.currslot == currslot);
		}
		inline bool
			operator!=(const iterator& x) const{
			return !operator == (x);
		}
	};
		/**************************************/
	//Node node1;
	//Node node2;
	//InnerNode innerNode_1;
	//InnerNode innerNode_2;
};