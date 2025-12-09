#include "emp-tool/emp-tool.h"
#include <iostream>
using namespace std;
using namespace emp;

void printt(block a) {
	//uint64_t i0 = _mm_extract_epi64(a, 0);
	//uint64_t i1 = _mm_extract_epi64(a, 1);
	//printf("%X %X\n", i0, i1);
	unsigned char *c = (unsigned char*)(&a);
	for(int i = 0; i < 16; ++i) printf("%b ", c[i]);
	printf("\n");
}





int main(void) {
	// sender
	block data[2], delta, table[2], w0, w1;
	MITCCRH<24> mi_gen;
	PRG prg;
	prg.random_block(&delta, 1);
	delta = delta | makeBlock(0x0, 0x1);
	delta = makeBlock(0x12345678ABCDABCD, 0x12345678ABCDABCD);
	// delta = makeBlock(0x010, 0x1); // just for testing
	mi_gen.setS(delta);
	halfblock deltaL, deltaR;
	splitBlock(delta, deltaL, deltaR);
	// receiver
	block data1[2];
	MITCCRH<24> mi_eva;
	mi_eva.setS(delta);
	block ret;
	bool z[5];
	halfblock halftable[3];

	cout << "Correctness ... ";

	prg.random_block(data, 2);
	// data[0] = data[1] = all_one_block;


	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
		{
			w0 = tregates_garble(data[0], data[1], delta,  halftable, z, &mi_gen, &prg);
			w1 = w0 ^ delta;
			// std::cout << w0 << std::endl << w1 << std::endl;
			if(i == 1) data1[0] = data[0] ^ delta; else data1[0] = data[0];
			if(j == 1) data1[1] = data[1] ^ delta; else data1[1] = data[1];
			ret = tregates_eval(data1[0], data1[1], halftable, z, &mi_eva);
			block ret1 = w0;
			if(i == 1 && j == 1) ret1 = w1;

			if(cmpBlock(&ret, &ret1, 1) == false)
			{
				cout << i << j << endl;
				cout << getLSB(data[0]) << getLSB(data[1]) << endl;
				cout << w0 << endl;
				cout << w1 << endl;
				cout << ret << endl;
				{cout << "wrong" << endl; abort();}
			}
			if(cmpBlock(&ret, &ret1, 1) == false) {std::cout << i << j << "wrong" << std::endl; abort();}
		}
	cout << "check\n";

	return 0;
}



