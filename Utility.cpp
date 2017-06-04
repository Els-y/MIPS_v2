#include <cstdio>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "Utility.h"

using namespace std;

bool read_data(
        int n,
        int d,
        float** &data,
        const char* file_name)
{
    FILE* fin = fopen(file_name, "r");
    if (!fin) {
        printf("%s doesn't exist!\n", file_name);
        return false;
    }

    int id;
    data = new float*[n];
    for (int i = 0; i < n; i++) {
        data[i] = new float[d];
        fscanf(fin, "%d", &id);
        for (int j = 0; j < d; j++) {
            fscanf(fin, "%f", &data[i][j]);
        }
    }

    printf("Finish reading %s\n", file_name);
    fclose(fin);

    return true;
}

bool read_data_with_id(int n, int d, float** & data, const char* file_name) {
    FILE* fin = fopen(file_name, "r");
    if (!fin) {
        printf("%s doesn't exist!\n", file_name);
        return false;
    }

    int id;
    data = new float*[n]; 
    for (int i = 0; i < n; i++) {
        data[i] = new float[d + 1]; // data[i][0] is dataID
        fscanf(fin, "%f", &data[i][0]);
        for (int j = 1; j <= d; j++) {
            fscanf(fin, "%f", &data[i][j]);
        }
    }

    printf("Finish reading %s\n", file_name);
    fclose(fin);

    return true;
}
/*
slot-info
1(int) number of slots(char)
*/

int getIndexSlotSize(int n, int d) {
    return 2 * (n + 1) * sizeof(int) +  (1 + d) * sizeof(float);
}

int getDataSlotSize(int d) {
    return 2 * sizeof(int) +  (1 + d + 20 * (d + 1)) * sizeof(float);
}
/*
index
dimension(d)  number(n)   rids        radius    center
1(int)        1(int)      n * 2(int)  1(float)  d(float)
*/
// n is n children
int calcSlotNumberOfIndex(int n, int d) {
    int slotSize = getIndexSlotSize(n, d);
    int slotInfoSize = 0, numberOfSlot = 0, total = 0;
    while (total < 64 * 1000) {
        numberOfSlot++;
        slotInfoSize = 1 * sizeof(int) + numberOfSlot * sizeof(char);
        total = slotSize * numberOfSlot + slotInfoSize;
    }
    // now, numberOfSlot is a little bigger than 64K, so return numberOfSlot-1
    return (numberOfSlot - 1);
}

/*
data
dimension(d)   number(n)  radius    center     data
1(int)         1(int)     1(float)  d(float)   n*(d+1)(float)
*/
// n is n data(vectors)
int calcSlotNumberOfData(int n, int d) {
    int slotSize = getDataSlotSize(d);
    int slotInfoSize = 0, numberOfSlot = 0, total = 0;
    while (total < 64 * 1000) {
        numberOfSlot++;
        slotInfoSize = 1 * sizeof(int) + numberOfSlot * sizeof(char);
        total = slotSize * numberOfSlot + slotInfoSize;
    }
    // now, numberOfSlot is a little bigger than 64K, so return numberOfSlot-1
    return (numberOfSlot - 1);
}
