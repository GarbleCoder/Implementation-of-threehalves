#include <iostream>

#include "emp-tool/emp-tool.h"
using namespace std;
using namespace emp;

class AbandonIO: public IOChannel<AbandonIO> { public:
	void send_data_internal(const void * data, int len) {
	}

	void recv_data_internal(void  * data, int len) {
	}
};

int port, party;
template <typename T>
void testHalfgateGarble(T* netio) {
	block* a = new block[128];
	block* b = new block[128];
	block* c = new block[128];

	PRG prg;
	prg.random_block(a, 128);
	prg.random_block(b, 128);

	string file = "./emp-tool/circuits/files/bristol_format/AES-non-expanded.txt";
	// string file = "./emp-tool/circuits/files/bristol_fashion/aes_256.txt";
	BristolFormat cf(file.c_str());

	if (party == BOB) {
		HalfGateEva<T>::circ_exec = new HalfGateEva<T>(netio);
		for (int i = 0; i < 100; ++i)
			cf.compute(c, a, b);
		delete HalfGateEva<T>::circ_exec;
	} else {
		AbandonIO* aio = new AbandonIO();
		HalfGateGen<AbandonIO>::circ_exec = new HalfGateGen<AbandonIO>(aio);

		auto start = clock_start();
		for (int i = 0; i < 100; ++i) {
			cf.compute(c, a, b);
		}
		double interval = time_from(start);
		cout << "Pure AES garbling speed : " << 100 * 6800 / interval << " million gate per second\n";
		delete aio;
		delete HalfGateGen<AbandonIO>::circ_exec;

		MemIO* mio = new MemIO();
		HalfGateGen<MemIO>::circ_exec = new HalfGateGen<MemIO>(mio);

		start = clock_start();
		for (int i = 0; i < 20; ++i) {
			mio->clear();
			for (int j = 0; j < 5; ++j)
				cf.compute(c, a, b);
		}
		interval = time_from(start);
		cout << "AES garbling + Writing to Memory : " << 100 * 6800 / interval << " million gate per second\n";
		delete mio;
		delete HalfGateGen<MemIO>::circ_exec;

		HalfGateGen<T>::circ_exec = new HalfGateGen<T>(netio);

		start = clock_start();
		for (int i = 0; i < 100; ++i) {
			cf.compute(c, a, b);
		}
		interval = time_from(start);
		cout << "AES garbling + Loopback Network : " << 100 * 6800 / interval << " million gate per second\n";

		delete HalfGateGen<T>::circ_exec;
	}

	delete[] a;
	delete[] b;
	delete[] c;
}




template <typename T>
void testThreegateGarble(T* netio) {
	block* a = new block[128];
	block* b = new block[128];
	block* c = new block[128];

	PRG prg;
	prg.random_block(a, 128);
	prg.random_block(b, 128);

	string file = "./emp-tool/circuits/files/bristol_format/AES-non-expanded.txt";
	BristolFormat cf(file.c_str());

	if (party == BOB) {
        TreGateEva<T>::circ_exec = new TreGateEva<T>(netio);
		// HalfGateEva<T>::circ_exec = new HalfGateEva<T>(netio);
		for (int i = 0; i < 100; ++i)
			cf.compute(c, a, b);
		delete TreGateEva<T>::circ_exec;
	} else {
		AbandonIO* aio = new AbandonIO();
        // HalfGateGen<AbandonIO>::circ_exec = new HalfGateGen<AbandonIO>(aio);
		TreGateGen<AbandonIO>::circ_exec = new TreGateGen<AbandonIO>(aio);
        
		auto start = clock_start();
		for (int i = 0; i < 100; ++i) {
			cf.compute(c, a, b);
		}
		double interval = time_from(start);
		cout << "Pure AES garbling speed : " << 100 * 6800 / interval << " million gate per second\n";
		delete aio;
		delete TreGateGen<AbandonIO>::circ_exec;

		MemIO* mio = new MemIO();
		TreGateGen<MemIO>::circ_exec = new TreGateGen<MemIO>(mio);

		start = clock_start();
		for (int i = 0; i < 20; ++i) {
			mio->clear();
			for (int j = 0; j < 5; ++j)
				cf.compute(c, a, b);
		}
		interval = time_from(start);
		cout << "AES garbling + Writing to Memory : " << 100 * 6800 / interval << " million gate per second\n";
		delete mio;
		delete TreGateGen<MemIO>::circ_exec;

		
		TreGateGen<T>::circ_exec = new TreGateGen<T>(netio);
		start = clock_start();
		for (int i = 0; i < 100; ++i) {
			cf.compute(c, a, b);
		}
		interval = time_from(start);
		cout << "AES garbling + Loopback Network : " << 100 * 6800 / interval  << " million gate per second\n";
		delete TreGateGen<T>::circ_exec;


	}

	delete[] a;
	delete[] b;
	delete[] c;
}













int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);


	// cout << "ThreeGates\n";
	// NetIO* ThreeNetio = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
	// testThreegateGarble<NetIO>(ThreeNetio);
	// delete ThreeNetio;

	cout << "ThreeGates\n";
	HighSpeedNetIO* ThreeNetio = new HighSpeedNetIO(party == ALICE ? nullptr : "127.0.0.1", port, port+1);
	testThreegateGarble<HighSpeedNetIO>(ThreeNetio);
	delete ThreeNetio;


	cout << "HalfGates\n";
	HighSpeedNetIO* HalfNetio = new HighSpeedNetIO(party == ALICE ? nullptr : "127.0.0.1", port, port+1);
	testHalfgateGarble<HighSpeedNetIO>(HalfNetio);
	delete HalfNetio;
}













































// #include <iostream>

// #include "emp-tool/emp-tool.h"
// using namespace std;
// using namespace emp;

// class AbandonIO: public IOChannel<AbandonIO> { public:
// 	void send_data_internal(const void * data, int len) {
// 	}

// 	void recv_data_internal(void  * data, int len) {
// 	}
// };

// int port, party;


// int main(int argc, char** argv) {
// 	parse_party_and_port(argv, &party, &port);
// 	cout << "Using NetIO\n";
// 	NetIO* netio = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
// 	test<NetIO>(netio);
// 	delete netio;

// 	// cout << "Using HighSpeedNetIO\n";
// 	// HighSpeedNetIO* hsnetio = new HighSpeedNetIO(party == ALICE ? nullptr : "127.0.0.1", port, port+1);
// 	// test<HighSpeedNetIO>(hsnetio);
// 	// delete hsnetio;
// }













