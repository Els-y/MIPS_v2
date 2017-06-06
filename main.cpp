#include <iostream>

#include "BallTree.h"
#include "Utility.h"

using std::cout;
using std::endl;

#define MNIST

#ifdef MNIST
char dataset[L] = "Mnist";
int n = 60000, d = 50;
int qn = 1000;

int delete_n = 100, insert_n = 100;

#endif

#ifdef NETFIX
char dataset[L] = "Netflix";
int n = 17770, d = 50;
int qn = 1000;
#endif

#ifdef test
char dataset[L] = "test";
int n = 2, d = 50;
int qn = 1;
#endif

int main() {
    char data_path[L], query_path[L];
    char index_path[L], output_path[L];
    float** data = nullptr;
    float** query = nullptr;

    sprintf(data_path, "%s/src/dataset.txt", dataset);
    sprintf(query_path, "%s/src/query.txt", dataset);
    sprintf(index_path, "%s/index", dataset);
    sprintf(output_path, "%s/dst/answer.txt", dataset);

    if (!read_data(n, d, data, data_path)) {
        return 1;
    }

    BallTree ball_tree1;
    ball_tree1.buildTree(n, d, data);
//    ball_tree1.buildQuadTree(n, d, data);
    ball_tree1.storeTree(index_path);

    if (!read_data(qn, d, query, query_path)) {
        return 1;
    }

    FILE* fout = fopen(output_path, "w");
    if (!fout) {
        printf("can't open %s!\n", output_path);
        return 1;
    }

    BallTree ball_tree2;
    ball_tree2.restoreTree(index_path);

    /*
        Test Delete
     */
//    char delete_path[L];
//    float** delete_data;
//
//    sprintf(delete_path, "%s/delete_src/deleteData.txt", dataset);
//    read_data_with_id(delete_n, d, delete_data, delete_path);
//
//    cout << "Deleting..." << endl;
//    for (int i = 0; i < delete_n; i++) {
//        ball_tree2.deleteData(d, delete_data[i]);
//    }
//    cout << "Delete Success!" << endl;
//
//    /*
//        Test Insert
//   //   */
//    char insert_path[L];
//    float** insert_data; // index 0 is ID
//
//    sprintf(insert_path, "%s/insert_src/insertData.txt", dataset);
//    read_data_with_id(insert_n, d, insert_data, insert_path);
//
//    cout << "Inserting..." << endl;
//    for (int i = 0; i < insert_n; i++) {
//        // cout << "insert " << insert_data[i][0] << endl;
//        ball_tree2.insertData(d, insert_data[i]);
//    }
//    cout << "Insert Success!" << endl;
    // -----------------------

    /*
        Query
     */

    cout << "MIP Searching..." << endl;
    for (int i = 0; i < qn; i++) {
        int index = ball_tree2.mipSearch(d, query[i]);
        fprintf(fout, "%d\n", index);
    }
    cout << "MIP Search Success!" << endl;
    fclose(fout);

    for (int i = 0; i < n; i++) {
        delete[] data[i];
    }

    for (int i = 0; i < qn; i++) {
        delete[] query[i];
    }

    
    /*
     * Test Query
     */
//    char queryTest_path[L], outputTest_path[L];
//    sprintf(queryTest_path, "%s/src/query.txt", dataset);
//    sprintf(outputTest_path, "%s/dst/testAnswer.txt", dataset);
//    float** queryTest = nullptr;
//    if(!read_data(qn, d, queryTest, queryTest_path));
//    FILE* test_fout = fopen(outputTest_path, "w");
//    if (!test_fout) {
//        printf("can't open %s!\n", outputTest_path);
//        return 1;
//    }
//    for (int j = 0; j < qn; ++j) {
//        int testIndex = ball_tree1.testQuery(n, d, data, queryTest[j]);
//        fprintf(test_fout, "%d\n", testIndex);
//    }
//    fclose(test_fout);
//
//    for (int k = 0; k < qn; ++k) {
//        delete[] queryTest[k];
//    }

    /*
        Test getFreeRid
     */
    // Rid freeRid = ball_tree2.getFreeRidForData(index_path);
    // cout << freeRid.page << ' ' << freeRid.slot << endl; 
    return 0;
}
