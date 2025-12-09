#ifndef EMP_TREGATE_EVA_H__
#define EMP_TREGATE_EVA_H__
#include "emp-tool/utils/utils.h"
#include "emp-tool/utils/mitccrh.h"
#include "emp-tool/execution/circuit_execution.h"
#include <iostream>
namespace emp {

// inline block tregates_eval(block A, block B, const halfblock *table, bool * z, MITCCRH<24> *mitccrh) {
// 	block W;
// 	bool sa, sb;
// 	halfblock halfInputLabels[4];
// 	splitBlock(A, halfInputLabels[0], halfInputLabels[1]);
// 	splitBlock(B, halfInputLabels[2], halfInputLabels[3]);
// 	sa = getLSB(A);
// 	sb = getLSB(B);

// 	block H[3];
// 	halfblock HashLabels[3];
// 	halfblock tmp;
// 	bool MSBH[3];
// 	H[0] = A;
// 	H[1] = B;
// 	H[2] = A ^ B;
// 	mitccrh->hash<3,1>(H);
// 	splitBlock(H[0], HashLabels[0], tmp);
// 	splitBlock(H[1], HashLabels[1], tmp);
// 	splitBlock(H[2], HashLabels[2], tmp);
// 	MSBH[0] = getMSB(H[0]); // obtain the first bit of the hash result to compute the coeffice.
// 	MSBH[1] = getMSB(H[1]);
// 	MSBH[2] = getMSB(H[2]);
// 	// std::cout << MSBH[0] << std::endl << MSBH[1] << std::endl << MSBH[2] << std::endl;
// 	//define V matrix
// 	// bool V00[10] = {1, 0, 0, 0, 0, // 2 * 5
// 	// 0, 1, 0, 0, 0};
// 	// bool V01[10] = {1, 0, 0, 0, 1,
// 	// 0, 1, 0, 1, 1};
// 	// bool V10[10] = {1, 0, 1, 0, 1,
// 	// 0, 1, 0, 0, 1};
// 	// bool V11[10] = {1, 0, 1, 0, 0,
// 	// 0, 1, 0, 1, 0};

// 	// bool* V_table[2][2] = {
//     // {V00, V01},
//     // {V10, V11}
// 	// };
// 	halfblock res1[2], res2[2], res3[2];
// 	// define Pr matrix
// 	bool Rp00[8] = {0, 0, 1, 0, 0, 1, 0, 0};// 2 * 4
// 	bool Rp01[8] = {0, 0, 1, 0, 0, 0, 0, 0};
// 	bool Rp10[8] = {0, 0, 0, 0, 0, 1, 0, 0};
// 	bool Rp11[8] = {0, 0, 0, 0, 0, 0, 0, 0};

// 	bool* Rp_table[2][2] = {
// 		{Rp00, Rp01},
// 		{Rp10, Rp11}
// 	};


// 	bool R[8];
// 	bool M[6] = {1, 0, 1, 0, 1, 1}; // 2 * 3
// 	// bool* V = V_table[sa][sb];
// 	halfblock G[5] = {0, 0, table[0], table[1], table[2]};
	
// 	bool r[2];
// 	bool tempbool1[2];
// 	bool tempbool2[2];
// 	// gf2MatMul_bool(tempbool1, V_table[sa][sb], z, 2, 5, 1);

// 	tempbool1[0] = z[0] ^ (sa & z[2]) ^ ((sa ^ sb) & z[4]);
// 	tempbool1[1] = z[1] ^ (sb & z[3]) ^ ((sa ^ sb) & z[4]);

// 	// gf2MatMul_bool(tempbool2, M, MSBH, 2, 3, 1);

// 	tempbool2[0] = MSBH[0] ^ MSBH[2];
// 	tempbool2[1] = MSBH[1] ^ MSBH[2];

// 	r[0] = tempbool1[0] ^ tempbool2[0];
// 	r[1] = tempbool1[1] ^ tempbool2[1];

// 	// std::cout << "z:" << z[0] << z[1] << z[2] << z[3] << z[4] << std::endl;
// 	// std::cout << "r:" << r[0] << r[1] << std::endl;

// 	// compute matrix R
// 	bool S1[8] = {1, 1, 1, 0, 1, 0, 0, 1};
// 	bool S2[8] = {1, 0, 0, 1, 0, 1, 1, 1};
// 	for(int i = 0; i < 8; i++)
// 		R[i] = (r[0] & S1[i]) ^ (r[1] & S2[i]) ^ Rp_table[sa][sb][i];
// 	// std::cout << "Ri" << std::endl;
// 	// for(int i = 0; i < 8; i++)
// 	// 	std::cout << R[i] << std::endl;
// 	// also need to xor Hash

// 	// gf2MatMul_bool_halfblock(res1, V_table[sa][sb], G, 2, 5, 1);

// 	// res1[0] = (sa & table[0]) ^ ((sa ^ sb) & table[2]);
// 	// res1[1] = (sb & table[1]) ^ ((sa ^ sb) & table[2]);	

// 	res1[0] = ((-(uint64_t) sa) & table[0]) ^ ((-(uint64_t) (sa ^ sb))& table[2]);
// 	res1[1] = ((-(uint64_t) sb) & table[1]) ^ ((-(uint64_t) (sa ^ sb))& table[2]);

// 	// std::cout << "res" << res1[0] << std::endl << res1[1] << std::endl;
// 	// gf2MatMul_bool_halfblock(res2, M, HashLabels, 2, 3, 1);
// 	res2[0] = HashLabels[0] ^ HashLabels[2];
// 	res2[1] = HashLabels[1] ^ HashLabels[2];
// 	gf2MatMul_bool_halfblock(res3 ,R, halfInputLabels, 2, 4, 1);
// 	W = makeBlock(res1[0]^res2[0]^res3[0], res1[1]^res2[1]^res3[1]);
// 	return W;
// }




inline block tregates_eval(block A, block B, const halfblock *table, bool * z, MITCCRH<24> *mitccrh) {
	block C;
	bool sa, sb;
	halfblock LA, RA, LB, RB, LC, RC;
	splitBlock(A, LA, RA);
	splitBlock(B, LB, RB);
	sa = getLSB(A);
	sb = getLSB(B);

	block H[3];
	halfblock HashLabels[3];
	bool MSBH[3];
	H[0] = A;
	H[1] = B;
	H[2] = A ^ B;
	// mitccrh->hash<3,1>(H);
	mitccrh->AesniTmmoBatch_FixedKey<3,1>(H);
	for(int i = 0; i < 3; i++)
	{
		HashLabels[i] = LBlock(H[i]);
		MSBH[i] = getMSB(H[i]);
	}

	/*
	V_{sa,sb} = 1, 0, sa, 0,  sa \oplus sb,
			    0, 1, 0,  sb, sa \oplus sb
	*/
	uint8_t Rp[4] = {0b00100100u, 0b00000100u, 0b00100000u, 0b00000000u};
	uint8_t S1 = 0b10010111u;
	uint8_t S2 = 0b11101001u;

	bool r1 = z[0] ^ (sa & z[2]) ^ ((sa ^ sb) & z[4]) ^ MSBH[0] ^ MSBH[2];
	bool r2 = z[1] ^ (sb & z[3]) ^ ((sa ^ sb) & z[4]) ^ MSBH[1] ^ MSBH[2];
	uint8_t uint8r1 = -r1;
	uint8_t uint8r2 = -r2;
	uint8_t R = (uint8r1 & S1) ^ (uint8r2 & S2) ^ Rp[2*sa+sb];
	LC = ((-(uint64_t) sa) & table[0]) ^ ((-(uint64_t) (sa ^ sb))& table[2]) ^ HashLabels[0] ^ HashLabels[2] ^ (-(R & uint8_t(0x1)) & LA) ^ (-((R >> 1) & uint8_t(0x1)) & RA) ^ (-((R >> 2) & uint8_t(0x1)) & LB) ^ (-((R >> 3) & uint8_t(0x1)) & RB);
	RC = ((-(uint64_t) sb) & table[1]) ^ ((-(uint64_t) (sa ^ sb))& table[2]) ^ HashLabels[1] ^ HashLabels[2] ^ (-((R >> 4)  & uint8_t(0x1)) & LA) ^ (-((R >> 5) & uint8_t(0x1)) & RA) ^ (-((R >> 6) & uint8_t(0x1)) & LB) ^ (-((R >> 7) & uint8_t(0x1)) & RB);


	// gf2MatMul_bool_halfblock(res3 ,R, halfInputLabels, 2, 4, 1);
	C = makeBlock(LC, RC);
	return C;
}









template<typename T>
class TreGateEva:public CircuitExecution {
public:
	T * io;
	block constant[2];
	MITCCRH<24> mitccrh;
	TreGateEva(T * io) :io(io) {
		set_delta();
		block tmp;
		io->recv_block(&tmp, 1);
		mitccrh.setS(tmp);
	}
	void set_delta() {
		io->recv_block(constant, 2);
	}
	block public_label(bool b) override {
		return constant[b];
	}

	// tregates_eval(block A, block B, const halfblock *table, bool * z, MITCCRH<8> *mitccrh) 

	block and_gate(const block& a, const block& b) override {
		halfblock table[3];
		bool z[8];
		io->recv_halfblock(table, 3);
		io->recv_bool(z, 5);
		return tregates_eval(a, b, table, z, &mitccrh);
	}
	block xor_gate(const block& a, const block& b) override {
		return a ^ b;
	}
	block not_gate(const block&a) override {
		return xor_gate(a, public_label(true));
	}
	uint64_t num_and() override {
		return mitccrh.gid/2;
	}
};
}
#endif// HALFGATE_EVA_H__
