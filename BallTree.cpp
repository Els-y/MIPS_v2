#include "BallTree.h"
#include "Utility.h"
#include <vector>
#include <list>
#include <cmath>
#include <queue>
#include <cmath>
#include <cstddef>
#include <string>
#include <cfloat>
#include <algorithm>
#include <cstring>
#include <iostream>
using namespace std;

BallTree::BallTree() {
    root = NULL;
}

BallTree::~BallTree() {
    clearTree(root);
}

void BallTree::clearTree(Node* root) {
    if (root != NULL) {
        vector<Node*> & childPtr = root->getPtr();
        for (auto p: childPtr) {
            clearTree(p);
        }
        delete root;
    }
}

int BallTree::getDimension() const {
    return dimension;
}

bool BallTree::buildTree(
        int n,
        int d,
        float** data) {

    list<vector<float> > listData = transform(n, d, data);

    dimension = 2;

    root = makeBallTree(listData);

    return true;
}

Node* BallTree::makeBallTree(list<vector<float> > & data) {
    int n = data.size();

    vector<float> center = getMean(data);
    vector<float> furthest = getFurthest(center, data);
    float radius = getDistance(center, furthest);

    Node* node = new Node();
    node->setRadius(radius);
    node->setCenter(center);

    if (n <= N0) {
        node->setType(1);
        list<vector<float> > & nodeData = node->getData();

        while (!data.empty()) {
            nodeData.push_back(data.front());
            data.pop_front();
        }

    } else {
        node->setType(0);
        vector<vector<float> > edgePoints = makeBallTreeSplit(data);
        vector<list<vector<float> > > subData = splitData(edgePoints, data);
        vector<list<vector<float> > > secondSubData;

        vector<Node*> & ptr = node->getPtr();
        if (dimension == 4) {
            for (int i = 0; i < 2; ++i) {
                vector<vector<float> > subEdgePoints = makeBallTreeSplit(subData[i]);
                vector<list<vector<float> > > tmpSubData = splitData(subEdgePoints, subData[i]);

                for (int j = 0; j < tmpSubData.size(); ++j) {
                    if (!tmpSubData[j].empty()) {
                        secondSubData.push_back(tmpSubData[j]);
                    }
                }
            }

            for (int i = 0; i < secondSubData.size(); ++i) {
                ptr.push_back(makeBallTree(secondSubData[i]));
            }
        } else {
            for (int i = 0; i < 2; ++i) {
                ptr.push_back(makeBallTree(subData[i]));
            }
        }
    }

    return node;
}

vector<vector<float> > BallTree::makeBallTreeSplit(const list<vector<float> > & data) {
    vector<float> pointX = data.front();
    vector<float> pointA = getFurthest(pointX, data);
    vector<float> pointB = getFurthest(pointA, data);
    vector<vector<float> > res(2);

    res[0] = pointA;
    res[1] = pointB;

    return res;
}

vector<list<vector<float> > > BallTree::splitData(const vector<vector<float> > & edgePoints, list<vector<float> > & data) {
    int size = edgePoints.size();
    vector<list<vector<float> > > subData(size, list<vector<float> >());

    while (!data.empty()) {
        vector<float> v = data.front();
        int index = findNeastIndex(v, edgePoints);

        subData[index].push_back(v);
        data.pop_front();
    }

    return subData;
}

// add data id and convert float** to list<vector<float> >
list<vector<float> > BallTree::transform(int n, int d, float** data) {
    list<vector<float> > res;

    for (int i = 0; i < n; ++i) {
        vector<float> v(1, i + 1);
        v.insert(v.end(), data[i], data[i] + d);
        res.push_back(v);
    }

    return res;
}

vector<float> BallTree::getFurthest(const vector<float> & pointX, const list<vector<float> > & data) {
    vector<float> res;

    float maxDist = -1;

    for (auto tmp: data) {
        float curDist = getDistance(pointX, tmp);

        if (curDist > maxDist) {
            maxDist = curDist;
            res = tmp;
        }
    }

    return res;
}

float BallTree::getDistance(const vector<float> & pointA, const vector<float> & pointB) {
    float sum = 0;
    int length = pointA.size();

    for (int i = 1; i < length; ++i) {
        sum += (pointA[i] - pointB[i]) * (pointA[i] - pointB[i]);
    }

    return sqrt(sum);
}

vector<float> BallTree::getMean(const list<vector<float> > & data) {
    int n = data.size();
    int size = data.front().size();

    vector<float> res(size, 0);
    res[0] = -1;

    for (auto v: data) {
        for (int i = 1; i < size; ++i) {
            res[i] += v[i];
        }
    }

    for (int i = 1; i < size; ++i) {
        res[i] /= n;
    }

    return res;
}

int BallTree::findNeastIndex(const vector<float> & pointX, const vector<vector<float> > & edgePoints) {
    int size = edgePoints.size();

    float minDist = -1;
    int index = -1;

    for (int i = 0; i < size; ++i) {
        float curDist = getDistance(pointX, edgePoints[i]);

        if (index == -1 || curDist < minDist) {
            minDist = curDist;
            index = i;
        }
    }

    return index;
}

// restore the whole tree except the data in the origin tree
bool BallTree::restoreTree(const char* index_path) {
    this->index_path = index_path;
    root = readNode(index_path, Rid(0, 0));
    dimension = root->getChildren().size();

    // restore tree without data
    queue<Node*> que;
    que.push(root);
    while (!que.empty()) {
        Node* tmp = que.front();
        que.pop();

        std::vector<Rid>& vRid = tmp->getChildren();
        for (auto rid : vRid) {
            if (rid.page < 1000) { // index node
                Node* tmpIndex = readNode(index_path, rid);
                tmp->getPtr().push_back(tmpIndex);
                que.push(tmpIndex);
            } else {
                tmp->getPtr().push_back(readDataNodeWithoutData(index_path, rid));
            }
        }
    }
    return root != nullptr;
}

Node* BallTree::readDataNodeWithoutData(const char * index_path, const Rid & rid) {
    char path[64];
    sprintf(path, "%s/%d", index_path, rid.page);
    ifstream in;
    in.open(path, ios::in|ios::binary);

    Node* retNode = new Node();
    retNode->setRid(rid);

    readData(in, rid.slot);

    int d, n, tmp;
    in.read((char *) &d, sizeof(int));
    in.read((char *) &n, sizeof(int));

    float tmpf;
    // radius
    in.read((char *) &tmpf, sizeof(float));
    retNode->setRadius(tmpf);

    // center
    retNode->getCenter().push_back((float)-1);
    for (int i = 0; i < d; i++) {
        in.read((char *) &tmpf, sizeof(float));
        retNode->getCenter().push_back(tmpf);
    }

    retNode->setType(1);
    in.close();
    return retNode;
}

float BallTree::getMIP(Query& q, Node* T) {
    return getInnerProduct(q.getData(), T->getCenter())
           + T->getRadius() * getNorm(q.getData());
}

float BallTree::getNorm(const vector<float> & q) {
    vector<float> zero(q.size(), 0);
    return getDistance(q, zero);
}

bool BallTree::isLeaf(Node* T) {
    return T->getChildren().empty();
}

void BallTree::linearSearch(Query& query, const list<vector<float> > & S) {
    for (auto each : S) {
        float dis = getInnerProduct(query.getData(), each);
        if(dis >= query.getMaxInnerProduct()) {
            query.setBM(int(each[0]));
            query.setMaxInnerProduct(dis);
        }
    }
}

float BallTree::getInnerProduct(const vector<float> & pointA,
                                const vector<float> & pointB) {
    float sum = 0;
    int length = (int)pointA.size();

    for (int i = 1; i < length; ++i) {
        sum += pointA[i] * pointB[i];
    }

    return sum;
}


void BallTree::treeSearch(Query& query, Node* root) {
    if (query.getMaxInnerProduct() < getMIP(query, root)) {
        if(isLeaf(root)) {
            this->linearSearch(query, root->getData());
        } else {
            // 保证vchildren中Rid的顺序和vptr中Node*的顺序一致，即对应同一个Node.
            // 利用vIndex记录得到排序后的vptr在原来root->getPtr()数组，
            // 也即children数组的位置
            std::vector<Rid>& vchildren = root->getChildren();
            std::vector<Node* > vptr = root->getPtr();
            sort(vptr.begin(), vptr.end(), [&](Node* a, Node* b) {
                return getMIP(query, a) > getMIP(query, b);
            });
            std::vector<int> vIndex;
            for (int i = 0; i < vptr.size(); i++) {
                for (int j = 0; j < vptr.size(); j++) {
                    if (vptr[i] == root->getPtr()[j]) {
                        vIndex.push_back(j);
                        break;
                    }
                }
            }
            for (int i = 0; i < vptr.size(); i++) {
                if (isLeaf(vptr[i])) {
                    Node* dataNode = readNode(this->index_path.c_str(), vchildren[vIndex[i]]);
                    treeSearch(query, dataNode);
                    delete dataNode;
                } else {
                    treeSearch(query, vptr[i]);
                }
            }
        }
    }
}

// query has no id!!
int BallTree::mipSearch(int d, float* query) {
    vector<float> v(1, -1);
    v.insert(v.end(), query, query + d);
    Query q(v);
    treeSearch(q, root);
    return q.getBM();
}

bool BallTree::buildQuadTree(
        int n,
        int d,
        float** data) {

    list<vector<float> > listData = transform(n, d, data);

    dimension = 4;

    root = makeBallTree(listData);

    return true;
}

bool BallTree::storeTree(const char* index_path) {
    this->index_path = index_path;

    // set rid for every node using BFS
    setRid(root);

    // clarify index node's children(in rid form)
    modifyRid(root);

    // store node in file
    storeInFile(index_path, root);

    return true;
}

void BallTree::setRid(Node* root) {
    // dimension of the vector
    int d = root->getCenter().size() - 1;

    // number of children
    int n = root->getPtr().size();

    // how many records one page can hold
    int numberOfIndex = calcSlotNumberOfIndex(n, d);
    int numberOfData = calcSlotNumberOfData(N0, d);

    int indexPid = 0, indexSid = 0, dataPid = 1000, dataSid = 0;

    queue<Node*> que;
    que.push(root);
    while (!que.empty()) {
        Node* tmp = que.front();
        que.pop();

        // judge it an index node or a data node
        if (tmp->isIndex()) { // an index node
            Rid rid(indexPid, indexSid);
            tmp->setRid(rid);
            if (indexSid == numberOfIndex - 1) {
                indexPid++;
                indexSid = 0;
            } else {
                indexSid++;
            }

            // find children
            std::vector<Node* > v = tmp->getPtr();
            for (int i = 0; i < v.size(); i++) {
                que.push(v[i]);
            }
        } else { // a data node
            Rid rid(dataPid, dataSid);
            tmp->setRid(rid);
            if (dataSid == numberOfData - 1) {
                dataPid++;
                dataSid = 0;
            } else {
                dataSid++;
            }

        }
    }
}

void BallTree::modifyRid(Node* root) {
    queue<Node*> que;
    que.push(root);
    while (!que.empty()) {
        Node* tmp = que.front();
        que.pop();

        // judge it an index node or a data node
        if (tmp->isIndex()) { // an index node
            // find children
            std::vector<Node* > v = tmp->getPtr();
            for (int i = 0; i < v.size(); i++) {
                tmp->getChildren().push_back(v[i]->getRid());
                que.push(v[i]);
            }
        }
    }
}

void BallTree::storeInFile(const char* index_path, Node* root) {
    // dimension of the vector
    int d = root->getCenter().size() - 1;

    // number of children
    int n = root->getPtr().size();

    // how many records one page can hold
    int numberOfIndex = calcSlotNumberOfIndex(n, d);
    int numberOfData = calcSlotNumberOfData(N0, d);

    int iPid = 0, iSid = 0, dPid = 1000, dSid = 0;
    createNewPage(index_path, iPid);
    createNewPage(index_path, dPid);
    queue<Node*> que;
    que.push(root);
    while (!que.empty()) {
        Node* tmp = que.front();
        que.pop();

        // judge it an index node or a data node
        if (tmp->isIndex()) { // an index node
            storeNode(index_path, *tmp);

            if (iSid == numberOfIndex - 1) {
                char path[64];
                sprintf(path, "%s/%d", index_path, iPid);
                writeSlotInfo(path, numberOfIndex, iSid + 1);

                iPid++;
                iSid = 0;
                createNewPage(index_path, iPid);
            } else {
                iSid++;
            }

            // find children
            std::vector<Node* > v = tmp->getPtr();
            for (int i = 0; i < v.size(); i++) {
                que.push(v[i]);
            }
        } else { // a data node
            storeNode(index_path, *tmp);

            if (dSid == numberOfData - 1) {
                char path[64];
                sprintf(path, "%s/%d", index_path, dPid);
                writeSlotInfo(path, numberOfData, dSid + 1);
                dPid++;
                dSid = 0;
                createNewPage(index_path, dPid);
            } else {
                dSid++;
            }
        }
    }

    // Fill the rest unused slots
    if (iSid != 0) {
        char path[64];
        sprintf(path, "%s/%d", index_path, iPid);

        appendInvalidSlots(path, getIndexSlotSize(n, d), numberOfIndex - iSid);

        writeSlotInfo(path, numberOfIndex, iSid);
    }

    if (dSid != 0) {
        char path[64];
        sprintf(path, "%s/%d", index_path, dPid);

        appendInvalidSlots(path, getDataSlotSize(d), numberOfData - dSid);

        writeSlotInfo(path, numberOfData, dSid);
    }
}

void BallTree::appendInvalidSlots(const char* path, int slotSize, int slotNumber) {
    ofstream out;
    out.open(path, ios::out | ios::binary | ios::app);

    char temp[slotSize - sizeof(int) * 2];

    memset(temp, '0', sizeof(temp));

    for (int i = 0; i < slotNumber; i++) {
        int d = 50;
        int tempn = 0;

        out.write((char*)&d, sizeof(int));
        out.write((char*)&tempn, sizeof(int));

        out.write(temp, sizeof(temp));
    }
    out.close();
}

void BallTree::storeDataNode(ofstream& out, Node& node) {
    int d = node.getCenter().size() - 1,
            n = node.getData().size();
    // dimension
    out.write((char* ) &d, sizeof(int));
    // n cha
    out.write((char* ) &n, sizeof(int));

    // radius
    float r = node.getRadius();
    out.write((char* ) &r, sizeof(float));

    // center
    std::vector<float> vcenter = node.getCenter();
    for (int i = 1; i < vcenter.size(); i++) {
        out.write((char* ) &vcenter[i], sizeof(float));
    }

    list<std::vector<float> > ldata = node.getData();
    // data
    for (int i = 0; i < n; i++) {
        std::vector<float> vf = ldata.front();
        ldata.pop_front();
        for (int i = 0; i < vf.size(); i++) {
            out.write((char* ) &vf[i], sizeof(float));
        }
    }

    // Fill the unused data place
    for (int i = 0; i < N0 - n; i++) {
        vector<float> invalidData(d + 1, 0);
        invalidData[0] = -1; // Invalid data ID
        for (int i = 0; i < invalidData.size(); i++) {
            out.write((char*) &invalidData[i], sizeof(float));
        }
    }
}

bool BallTree::storeNode(const char* index_path, Node & node) {
    char path[64];
    Rid rid = node.getRid();
    sprintf(path, "%s/%d", index_path, rid.page);

    ofstream out;
    out.open(path, ios::out|ios::binary|ios::app);

    if (node.isIndex()) {  // an index node
        storeIndexNode(out, node);
    } else {  //  a data node
        storeDataNode(out, node);
    }
    out.close();
    return true;
}

bool BallTree::writeSlotInfo(const char * path, int slotSize, int occupySize) {
    ofstream out;
    out.open(path, ios::out|ios::binary|ios::app);

    out.write((char* ) &slotSize, sizeof(int));

    char ch = '1';
    for (int i = 0; i < occupySize; i++) {
        out.write(&ch, sizeof(char));
    }
    char ch2 = '0';
    for (int i = 0; i < slotSize - occupySize; i++) {
        out.write(&ch2, sizeof(char));
    }
    out.close();
    return true;
}

void BallTree::createNewPage(const char* index_path, int pid) {
    char path[64];
    sprintf(path, "%s/%d", index_path, pid);
    ofstream out;
    out.open(path, ios::out|ios::binary);
    out.close();
}

Node* BallTree::readNode(const char * index_path, const Rid & rid) {
    char path[64];
    sprintf(path, "%s/%d", index_path, rid.page);
    ifstream in;
    in.open(path, ios::in|ios::binary);

    Node* retNode = new Node();
    retNode->setRid(rid);

    if (rid.page < 1000) {  // Index Node
        readIndex(in, rid.slot);

        int d, n, tmp1, tmp2;
        in.read((char *) &d, sizeof(int));
        in.read((char *) &n, sizeof(int));

        // children
        for (int i = 0; i < n; i++) {
            in.read((char *) &tmp1, sizeof(int));
            in.read((char *) &tmp2, sizeof(int));
            Rid child(tmp1, tmp2);
            retNode->getChildren().push_back(child);
        }

        for (int i = 0; i < dimension - n; i++) {
            in.seekg((int)in.tellg() + sizeof(int) * 2);
        }

        float tmpf;
        // radius
        in.read((char *) &tmpf, sizeof(float));
        retNode->setRadius(tmpf);

        // center
        retNode->getCenter().push_back((float)-1);
        for (int i = 0; i < d; i++) {
            in.read((char *) &tmpf, sizeof(float));
            retNode->getCenter().push_back(tmpf);
        }
        retNode->setType(0);
    } else {    // Data Node
        readData(in, rid.slot);

        int d, n, tmp;
        in.read((char *) &d, sizeof(int));
        in.read((char *) &n, sizeof(int));

        float tmpf;
        // radius
        in.read((char *) &tmpf, sizeof(float));
        retNode->setRadius(tmpf);

        // center
        retNode->getCenter().push_back((float)-1);
        for (int i = 0; i < d; i++) {
            in.read((char *) &tmpf, sizeof(float));
            retNode->getCenter().push_back(tmpf);
        }

        // data
        for (int i = 0; i < n; i++) {
            std::vector<float> v;
            for (int i = 0; i < d + 1; i++) {
                in.read((char *) &tmpf, sizeof(float));
                v.push_back(tmpf);
            }
            retNode->getData().push_back(v);
        }
        retNode->setType(1);
    }

    in.close();
    return retNode;
}

void BallTree::readIndex(ifstream& in, int num) {
    int d;
    in.read((char*) &d, sizeof(int));

    int totalSize = sizeof(int) * (2 + dimension * 2) + sizeof(float) * (1 + d);
    in.seekg((int)in.tellg() + totalSize * num - sizeof(int));
}

// Skip the size of Data in ifstream
void BallTree::readData(ifstream & in, int n) {
    int d;
    in.read((char* ) &d, sizeof(int));
    int totalSize = sizeof(int) * 2 + sizeof(float) * (1 + d + N0 * (d + 1));
    in.seekg((int)in.tellg() + totalSize * n - sizeof(int));
}

Rid BallTree::getFreeRid(const char * index_path, int pid) {
    char path[64];
    int d, n;
    while (true) {
        sprintf(path, "%s/%d", index_path, pid);

        string slotValidMap = getSlotValidMap(index_path, pid);

        if (slotValidMap == "NULL") {
            return Rid(pid, -1);
        }

        int sid = 0, result = -1;
        while (sid < slotValidMap.size()) {
            if (slotValidMap[sid] == '0') {
                result = sid;
                return Rid(pid, sid);
            } else {
                sid++;
            }
        }

        if (result == -1) {
            pid++;
        }
    }
}

string BallTree::getSlotValidMap(const char* index_path, int pid) {
    char path[64];
    int d, n;
    int slotNumber, slotSize;

    sprintf(path, "%s/%d", index_path, pid);

    ifstream in;
    in.open(path, ios::in | ios::binary);

    if (in.fail()) {
        return "NULL";
    }
    in.read((char*) &d, sizeof(int));
    in.read((char*) &n, sizeof(int));

    if (pid < 1000) {
        slotNumber = calcSlotNumberOfIndex(n, d);
        slotSize = getIndexSlotSize(n, d);
    } else {
        slotNumber = calcSlotNumberOfData(n, d);
        slotSize = getDataSlotSize(d);
    }

    in.seekg(0, in.end);
    int pageLength = in.tellg();

    int mapLength = pageLength - slotSize * slotNumber;

    char* map = new char[mapLength];

    in.seekg(slotNumber * slotSize);

    int slotNum;
    in.read((char*)& slotNum, sizeof(int));  // extract the first 4 bytes
    in.read(map, mapLength - sizeof(int));

    string mapStr = "";
    for (int i = 0; i < mapLength - sizeof(int); ++i) {
        mapStr += map[i];
    }

    in.close();
    return mapStr;
}

Rid BallTree::getFreeRidForIndex(const char * index_path) {
    return getFreeRid(index_path, 0);
}


Rid BallTree::getFreeRidForData(const char * index_path) {
    return getFreeRid(index_path, 1000);
}

void BallTree::createNewPage(const char* index_path, int pid, int slotNum, int slotSize) {
    cout << "create new page " << pid << ' ' << slotNum << ' ' << slotSize << endl;

    char path[64];
    sprintf(path, "%s/%d", index_path, pid);
    ofstream out;
    out.open(path, ios::out|ios::binary);
    out.close();

    appendInvalidSlots(path, slotSize, slotNum);
    writeSlotInfo(path, slotNum, 0);
}

// optional

// If the distance between data and node's center is bigger than node's radius,
// set its radius to this distance
void BallTree::updateNodeRadius(Node* node, vector<float> data) {
    float radius = node->getRadius();
    vector<float> center = node->getCenter();
    float dist = getDistance(center, data);

    if (dist > radius) {
        node->setRadius(dist);
    }
}


bool BallTree::insertData(int d, float* data) {
    vector<float> vec_data(data, data + d + 1);
    insertDataToIndexNode(root, vec_data);
    return true;
}

// Recursively insert the data to index node
void BallTree::insertDataToIndexNode(Node* indexNode, vector<float>& data) {
    updateNodeRadius(indexNode, data);

    vector<Rid> children = indexNode->getChildren();

    float minDist = FLT_MAX;
    Node* targetChild;

    for (int i = 0; i < children.size(); i++) {
        Node* childNode = readNode(this->index_path.c_str(), children[i]);

        float dist = getDistance(childNode->getCenter(), data);

        // if (dist < childNode->getRadius()) {
        //     targetChild = childNode;
        //     break;
        // }

        if (dist < minDist) {
            minDist = dist;
            targetChild = childNode;
        } else {
            delete childNode;
        }
    }

    if (targetChild->isIndex()) {
        insertDataToIndexNode(targetChild, data);
    } else {
        insertDataToDataNode(indexNode, targetChild, data);
    }

    delete targetChild;
}

// Calculate the total number of data in this node
// If the node is an index, the result is the sum of the data of its children.
int BallTree::getDataNum(Node* node) {
    if (node->isData()) {
        return node->getData().size();
    }

    // Index
    int totalNum = 0;
    vector<Rid> children = node->getChildren();
    for (int i = 0; i < children.size(); i++) {
        Node* tempNode = readNode(this->index_path.c_str(), children[i]);
        int temp = getDataNum(tempNode);

        totalNum += temp;
        delete tempNode;
    }
    return totalNum;
}

// Add the data to data node.
void BallTree::insertDataToDataNode(Node* fatherNode, Node* dataNode, vector<float>& data) {
    updateNodeRadius(dataNode, data);

    list<vector<float> > & dataList = dataNode->getData();

    if (dataList.size() < N0) {
        dataList.push_back(data);
        updateSlot(true, dataNode);
        return;
    }

    // no valid space
    splitDataNode(fatherNode, dataNode, data);
}

// Split this data node to an new index node
void BallTree::splitDataNode(Node* fatherNode, Node* dataNode, vector<float>& data) {
    list<vector<float> > newDataList = dataNode->getData();

    newDataList.push_back(data);
    Node* newIndex = makeBallTree(newDataList);

    setNewIndexNodesRid(newIndex);

    vector<Rid> & fatherChildren = fatherNode->getChildren();

    for (int i = 0; i < fatherChildren.size(); i++) {
        if (fatherChildren[i].slot == dataNode->getRid().slot && fatherChildren[i].page == dataNode->getRid().page) {
            fatherChildren[i] = newIndex->getRid();
            break;
        }
    }
    // delete old current data node
    updateSlot(false, dataNode);

    // update father node
    updateSlot(true, fatherNode);

    // store new index node
    updateSlot(true, newIndex);

    // store new data node
    vector<Node*> ptr = newIndex->getPtr();
    for (int i = 0; i < ptr.size(); i++) {
        updateSlot(true, ptr[i]);
    }

    delete newIndex;
}

// Set rid for the index node and its children
void BallTree::setNewIndexNodesRid(Node* node) {
    Rid newIndexRid = getFreeRidForIndex(this->index_path.c_str());

    while (newIndexRid.slot == -1) {
        createNewPage(this->index_path.c_str(), newIndexRid.page, getSlotNumber(node), getSlotSize(node));
        newIndexRid = getFreeRidForIndex(this->index_path.c_str());
    }

    node->setRid(newIndexRid);
    updateSlotInfo(true, node);

    vector<Node*> & ptr = node->getPtr();
    vector<Rid> & children = node->getChildren();

    for (int i = 0; i < ptr.size(); i++) {
        Rid newDataRid = getFreeRidForData(this->index_path.c_str());
        while (newDataRid.slot == -1) {
            createNewPage(this->index_path.c_str(), newDataRid.page, getSlotNumber(ptr[i]), getSlotSize(ptr[i]));
            newDataRid = getFreeRidForData(this->index_path.c_str());
        }

        ptr[i]->setRid(newDataRid);
        children.push_back(newDataRid);

        updateSlotInfo(true, ptr[i]);
    }
}

// valid == false, only make the slot map invalid
// valid == true, rewrite the slot
void BallTree::updateSlot(bool valid, Node* node) {
    Rid rid = node->getRid();

    char path[64];
    sprintf(path, "%s/%d", this->index_path.c_str(), rid.page);

    int slotSize = getSlotSize(node);

    ofstream out;
    out.open(path, ios::in | ios::out | ios::binary);

    out.seekp(slotSize * rid.slot);

    if (valid) {
        if (node->isData()) {
            // Data Node
            storeDataNode(out, *node);
        } else {
            // Index Node
            storeIndexNode(out, *node);
        }
    }
    out.close();
    updateSlotInfo(valid, node);
}

// Store the index node to its page and slot
void BallTree::storeIndexNode(ofstream& out, Node& node) {
    int d = node.getCenter().size() - 1,
            n = node.getChildren().size();

    out.write((char* ) &d, sizeof(int));
    // n cha
    out.write((char* ) &n, sizeof(int));

    // 2 * n rids
    std::vector<Rid> v = node.getChildren();
    for (int i = 0; i < n; i++) {
        out.write((char* ) &(v[i].page), sizeof(int));
        out.write((char* ) &(v[i].slot), sizeof(int));
    }

    Rid invalidRid(-1, -1);
    for (int i = 0; i < dimension - n; i++) {
        out.write((char*) &invalidRid.page, sizeof(int));
        out.write((char*) &invalidRid.slot, sizeof(int));
    }

    // radius
    float r = node.getRadius();
    out.write((char* ) &r, sizeof(float));

    // center
    std::vector<float> vcenter = node.getCenter();
    for (int i = 1; i < vcenter.size(); i++) {
        out.write((char* ) &vcenter[i], sizeof(float));
    }
}

// valid == false, make the slot valid '0'
// valid == true, make the slog invalid '1'
void BallTree::updateSlotInfo(bool valid, Node* node) {
    Rid rid = node->getRid();
    char path[64];
    sprintf(path, "%s/%d", this->index_path.c_str(), rid.page);

    ofstream out;
    out.open(path, ios::in | ios::out | ios::binary);
    out.seekp(getSlotSize(node) * getSlotNumber(node) + sizeof(int) + (rid.slot));

    char ch = valid ? '1' : '0';
    out.write(&ch, sizeof(char));
    out.close();
}

// Get the size of the slot when this node is stored in page
int BallTree::getSlotSize(Node* node) {
    int n, d;
    int slotSize, slotNumber;

    d = node->getCenter().size() - 1;
    if (node->isData()) {
        // Data Node
        return getDataSlotSize(d);
    } else {
        // Index Node
        n = dimension;
        return getIndexSlotSize(n, d);
    }
}

int BallTree::getSlotNumber(Node* node) {
    int n, d;
    int slotSize, slotNumber;

    d = node->getCenter().size() - 1;

    if (node->isData()) {
        // Data Node
        n = node->getData().size();
        return calcSlotNumberOfData(n, d);
    } else {
        // Index Node
        n = dimension;
        return calcSlotNumberOfIndex(n, d);
    }
}

bool BallTree::deleteData(int d, float* data) {
    vector<float> vec_data(data, data + d + 1);  // data has id
    deleteDataFromIndexNode(NULL, root, vec_data);
    return true;
}

// Recursively find data in this index node.
// If it's in, delete the data.
bool BallTree::deleteDataFromIndexNode(Node* father, Node* node, vector<float>& data) {
    bool flag = false;
    vector<Rid> children = node->getChildren();

    float minDist = FLT_MAX;
    Node* minChild = NULL;

    for (int i = 0; i < children.size(); i++) {
        if (flag) return true;
        Node* childNode = readNode(this->index_path.c_str(), children[i]);

        if (childNode->isIndex()) {
            flag = deleteDataFromIndexNode(node, childNode, data);
        } else {
            flag = deleteDataFromDataNode(father, node, childNode, data);
        }
        delete childNode;
    }
    delete minChild;
    return false;
}

// Find data in this data node.
// If it's in, delete the data.
bool BallTree::deleteDataFromDataNode(Node* grandfather, Node* father, Node* node, vector<float>& data) {
    list<vector<float> > & dataList = node->getData();

    for (auto i = dataList.begin(); i != dataList.end(); i++) {
        if (i->at(0) == data[0]) {
            dataList.erase(i);
            updateSlot(true, node);

            if (getDataNum(father) <= N0) {
                mergeIndexNode(grandfather, father);
            }
            cout << "delete data " << data[0] << " success" << endl;
            return true;
        }
    }
    return false;
}

// Merge this index node to a new data node
void BallTree::mergeIndexNode(Node* father, Node* node) {
    // delete this index node in index page
    updateSlot(false, node);

    vector<Rid> children = node->getChildren();
    list<vector<float> > totalData;

    for (int i = 0; i < children.size(); i++) {
        // delete the data node children
        Node* child = readNode(this->index_path.c_str(), children[i]);
        list<vector<float> > childData = child->getData();

        for (auto j = childData.begin(); j != childData.end(); j++) {
            totalData.push_back(*j);
        }
        // delete the data children in data page
        updateSlot(false, child);
        delete child;
    }

    // convert to data node
    node->getData() = totalData;
    node->getChildren().clear();
    node->getPtr().clear();
    node->setType(1);

    Rid oldRid = node->getRid();

    // store the converted data node to data page
    node->setRid(getFreeRidForData(this->index_path.c_str()));

    updateSlot(true, node);

    // update father's child
    vector<Rid> & fatherChildren = father->getChildren();
    for (int i = 0; i < fatherChildren.size(); i++) {
        if (fatherChildren[i].page == oldRid.page && fatherChildren[i].slot == oldRid.slot) {
            fatherChildren[i] = node->getRid();
        }
    }
    updateSlot(true, father);
}

// Node
void Node::setRid(const Rid & r) {
    rid = r;
}

Rid Node::getRid() const {
    return rid;
}

float Node::getRadius() const {
    return radius;
}

void Node::setRadius(const float & r) {
    radius = r;
}

vector<float> & Node::getCenter() {
    return center;
}

void Node::setCenter(const vector<float> & v) {
    center = v;
}

list<vector<float> > & Node::getData() {
    return data;
}

vector<Rid> & Node::getChildren() {
    return child;
}

vector<Node *> & Node::getPtr() {
    return ptr;
}

bool Node::isData() {
    return type == 1;
}

bool Node::isIndex() {
    return type == 0;
}

void Node::setType(int t) {
    type = t;
}
