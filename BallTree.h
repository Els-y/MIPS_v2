#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#define N0 20

#include <vector>
#include <list>
#include <fstream>

using std::vector;
using std::list;
using std::ifstream;
using std::ofstream;
using std::fstream;
using std::string;

// Node
struct Rid {
    int page;
    int slot;
    Rid(int p = 0, int s = 0) {
        page = p; slot = s;
    }
};

class Node {
private:
    Rid rid;
    float radius;
    vector<float> center;
    list<vector<float> > data;

    vector<Rid> child;
    vector<Node*> ptr;

    int type; // 1 for data, 0 for index

    void removeFromFile(string index_path);

public:
    void setRid(const Rid & r);

    bool isData();

    bool isIndex();

    void setType(int t);

    int getType() {return type;}

    Rid getRid() const;

    float getRadius() const;

    void setRadius(const float & r);

    vector<float> & getCenter();

    void setCenter(const vector<float> & v);

    list<vector<float> > & getData();

    vector<Rid> & getChildren();

    vector<Node *> & getPtr();
};

class Query {
private:
    vector<float> data;
    int bm; // ID of the best neighbor.
    float maxInnerProduct; // maximum inner product.
public:
    Query(const vector<float> & data) {
        this->data = data;
        bm = 0;
        maxInnerProduct = 0;
    }

    int getBM() {
        return bm;
    }

    void setBM(int bm) {
        this->bm = bm;
    }

    float getMaxInnerProduct() {
        return maxInnerProduct;
    }

    void setMaxInnerProduct(float maxInnerProduct) {
        this->maxInnerProduct = maxInnerProduct;
    }

    vector<float> getData() {
        return data;
    }
};

class BallTree {
private:
    int dimension;
    Node* root;
    string index_path;

    // buildTree
    Node* makeBallTree(list<vector<float> > & data);

    vector<vector<float> > makeBallTreeSplit(const list<vector<float> > & data);

    vector<list<vector<float> > > splitData(const vector<vector<float> > & edgePoints, list<vector<float> > & data);

    list<vector<float> > transform(int n, int d, float** data);

    vector<float> getFurthest(const vector<float> & pointX, const list<vector<float> > & data);

    float getDistance(const vector<float> & pointA, const vector<float> & pointB);

    vector<float> getMean(const list<vector<float> > & data);

    int findNeastIndex(const vector<float> & pointX, const vector<vector<float> > & edgePoints);

    void clearTree(Node* root);

    // store and read a tree
    Node* readDataNodeWithoutData(const char * index_path, const Rid & rid);

    Node* readNode(const char * index_path, const Rid & rid);

    bool storeNode(const char * index_path, Node & node);

    bool writeSlotInfo(const char * path, int slotSize, int occupySize);

    void storeInFile(const char* index_path, Node* root);

    void modifyRid(Node* root);

    void setRid(Node* root);

    void createNewPage(const char* index_path, int pid, int slotNum, int slotSize);

    void createNewPage(const char* index_path, int pid);

    void readIndex(ifstream & in, int n);

    // void readIndex(ifstream & in);

    void readData(ifstream & in, int n);

    void storeIndexNode(const char* path, int offset, Node* node);

    void storeDataNode(ofstream& out, Node& node);

    void storeIndexNode(ofstream& out, Node& node);

    void appendInvalidSlots(const char* path, int slotSize, int slotNumber);

    Rid getFreeRid(const char* index_path, int pid);

    string getSlotValidMap(const char* index_path, int pid);

    // Insert & Delete
    void updateNodeRadius(Node* node, vector<float> data);

    void insertDataToIndexNode(Node* indexNode, vector<float>& data);

    void insertDataToDataNode(Node* fatherNode, Node* dataNode, vector<float>& data);

    void splitDataNode(Node* fatherNode, Node* dataNode, vector<float>& data);

    int getSlotSize(Node* node);

    int getSlotNumber(Node* node);

    int getDataNum(Node* node);

    void updateSlot(bool valid, Node* node);

    void updateSlotInfo(bool valid, Node* node);

    void setNewIndexNodesRid(Node* node);

    bool deleteDataFromIndexNode(Node* father, Node* node, vector<float>& data);

    bool deleteDataFromDataNode(Node* grandfather, Node* father, Node* node, vector<float>& data);

    void mergeIndexNode(Node* father, Node* node);

    // Search
    float getMIP(Query& query, Node* T);

    float getNorm(const vector<float> & q);

    bool isLeaf(Node* T);

    void linearSearch(Query & q, const list<vector<float> > & S);

    void treeSearch(Query& query, Node* root);

    float getInnerProduct(const vector<float> & pointA,
                          const vector<float> & pointB);

public:
    BallTree();
    ~BallTree();

    Rid getFreeRidForIndex(const char * index_path);

    Rid getFreeRidForData(const char * index_path);

    int getDimension() const;

    bool buildTree(
            int n,
            int d,
            float** data);

    bool storeTree(
            const char* index_path);

    bool restoreTree(
            const char* index_path);

    int mipSearch(
            int d,
            float* query);

    // optional
    bool insertData(
            int d,
            float* data);

    // optional
    bool deleteData(
            int d,
            float* data);

    // optional
    bool buildQuadTree(
            int n,
            int d,
            float** data);
};

#endif



