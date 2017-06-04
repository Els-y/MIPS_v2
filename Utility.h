#ifndef __UTILITY_H
#define __UTILITY_H

#define L 256

bool read_data(
		int n,
		int d,
		float** &data,
		const char* file_name);

// preserve the id
bool read_data_with_id(int n, int d, float** & data, const char* file_name);

int getIndexSlotSize(int n, int d);
int getDataSlotSize(int d);
int calcSlotNumberOfIndex(int n, int d);
int calcSlotNumberOfData(int n, int d);

#endif


