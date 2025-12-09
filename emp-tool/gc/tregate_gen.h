#ifndef EMP_TREGATE_GEN_H__
#define EMP_TREGATE_GEN_H__
#include "emp-tool/utils/utils.h"
#include "emp-tool/utils/mitccrh.h"
#include "emp-tool/utils/prg.h"
#include "emp-tool/execution/circuit_execution.h"
#include "emp-tool/circuits/bit.h"
#include <iostream>
namespace emp {

/*
 * The half-gate garbling scheme, with improved hashing
 * [REF] Implementation of "Two Halves Make a Whole"
 * https://eprint.iacr.org/2014/756.pdf
 */

// 打印 64 位 halfblock
inline void printHalfblock(halfblock x, const std::string& name = "") {
    if (!name.empty()) std::cout << name << " = ";
    std::cout << std::bitset<64>(x) << std::endl;
}

// 打印 128 位 block
inline void printBlock(block b, const std::string& name = "") {
    alignas(16) uint64_t tmp[2];
    _mm_storeu_si128((__m128i*)tmp, b); // 拆成两个 64 位
    if (!name.empty()) std::cout << name << " = ";
    // 注意：tmp[1] 是高 64 位，tmp[0] 是低 64 位（小端）
    std::cout << std::bitset<64>(tmp[1]) << std::bitset<64>(tmp[0]) << std::endl;
}


inline block tregates_garble(block A0, block B0, block delta, halfblock *table, bool *z, MITCCRH<24> *mitccrh, PRG *prg) {
    bool pa = getLSB(A0);
    bool pb = getLSB(B0);  
    bool barpa = pa ^ 1;   //2 * barpa + \barpb is the position of 1 in the garbled table, in fact, bara and barb is the "a" and "b" in the paper "three halves makes a whole"
	bool barpb = pb ^ 1;

    //choose the coefficient of generation R_$
    bool rand[2]; 
    prg->random_bool(rand, 2); 

    uint32_t uint32rand0 = -uint32_t(rand[0]);
    uint32_t uint32rand1 = -uint32_t(rand[1]);


    static constexpr uint32_t VinvRp  = 0b00000000010010100100000010000100u; // Last 30 bit are the matrix V^-1 \times Rp
    static constexpr uint32_t VinvRa  = 0b00111110101001110111000000000000u; // 0b00111110101001110111000000000000
    static constexpr uint32_t VinvRb  = 0b00010111111110011001000000000000u; // 0b00010111111110011001000000000000
    static constexpr uint32_t VinvRr1 = 0b00010000110000100000001001000111u; // 0b00010000110000100000001001000111
    static constexpr uint32_t VinvRr2 = 0b00100000010000110000001110001001u; // 0b00100000010000110000001110001001
    static constexpr uint32_t VinvT[4] = {
    0b00000000110000110000100000010000u, // VinvT when pa = pb =0
    0b00000000110000000000000000000000u, // VinvT01 when pa = 0, pb = 1
    0b00010000000000110000000000000000u, // VinvT10 when pa = 1, pb = 0
    0b00010000000100000000000000000000u  // VinvT11 when pa = pb = 1
    };

    uint32_t uint32barpa = -uint32_t(barpa);
    uint32_t uint32barpb = -uint32_t(barpb);
    uint32_t VinvR = VinvRp ^ (uint32barpa & VinvRa) ^ (uint32barpb & VinvRb) ^ (uint32rand0 & VinvRr1) ^ (uint32rand1 & VinvRr2) ^ VinvT[static_cast<int>(2 * barpa + barpb)];


    block H[6];
    halfblock HashLabels[6];
    bool HF[6] = {};

    // Compute the hash function
	H[pa] = A0;
	H[barpa] = A0 ^ delta;
	H[2+pb] = B0;
	H[2+barpb] = B0 ^ delta;
    H[4+(pa^pb)] = A0 ^ B0;
    H[4+(pa^barpb)] = A0 ^ B0 ^ delta;

    halfblock LA, RA, LB, RB, LC, RC, LDelta, RDelta;
    splitBlock(H[0], LA, RA);
    splitBlock(H[2], LB, RB);
    splitBlock(delta, LDelta, RDelta); // split the input labels and the public offset


	// mitccrh->hash<3,2>(H); 
    mitccrh->AesniTmmoBatch_FixedKey<3,2>(H);
    for(int i = 0; i < 6; i++)
    {
        HF[i] = getMSB(H[i]);
        HashLabels[i] = LBlock(H[i]);
    }

    // Compute the outout label
    LC = (-(VinvR & uint64_t(0x1)) & LA) ^ (-((VinvR >> 1) & uint64_t(0x1)) & RA) ^ (-((VinvR >> 2) & uint64_t(0x1)) & LB) ^ (-((VinvR >> 3) & uint64_t(0x1)) & RB) ^ (-((VinvR >> 4) & uint64_t(0x1)) & LDelta) ^ (-((VinvR >> 5) & uint64_t(0x1)) & RDelta) ^ HashLabels[0] ^ HashLabels[4];
    VinvR = VinvR >> 6;
    RC = (-(VinvR & uint64_t(0x1)) & LA) ^ (-((VinvR >> 1) & uint64_t(0x1)) & RA) ^ (-((VinvR >> 2) & uint64_t(0x1)) & LB) ^ (-((VinvR >> 3) & uint64_t(0x1)) & RB) ^ (-((VinvR >> 4) & uint64_t(0x1)) & LDelta) ^ (-((VinvR >> 5) & uint64_t(0x1)) & RDelta) ^ HashLabels[2] ^ HashLabels[4];
    VinvR = VinvR >> 6;
    block C = makeBlock(LC, RC);
    // Compute the ciphertexts
    for(int i = 0; i < 3; i++)
    {
        table[i] = -(VinvR & uint64_t(0x1)) & LA; VinvR >>= 1;
        table[i] ^= -(VinvR & uint64_t(0x1)) & RA; VinvR >>= 1;
        table[i] ^= -(VinvR & uint64_t(0x1)) & LB; VinvR >>= 1;
        table[i] ^= -(VinvR & uint64_t(0x1)) & RB; VinvR >>= 1;
        table[i] ^= -(VinvR & uint64_t(0x1)) & LDelta; VinvR >>= 1;
        table[i] ^= -(VinvR & uint64_t(0x1)) & RDelta; VinvR >>= 1;
        table[i] ^= HashLabels[2 * i] ^ HashLabels[2 * i + 1];
    }

    // Computet the Dicing bits
    z[0] = rand[0] ^ HF[0] ^ HF[4];
    z[1] = rand[1] ^ HF[2] ^ HF[4];
    z[2] = barpa ^ HF[0] ^ HF[1];
    z[3] = barpb ^ HF[2] ^ HF[3];
    z[4] = barpa ^ barpb ^ HF[4] ^ HF[5];
	return C;
}


// inline block tregates_garble(block A0, block B0, block delta, halfblock *table, bool *z, MITCCRH<24> *mitccrh, PRG *prg) {
//     bool pa = getLSB(A0);
// 	bool pb = getLSB(B0);

//     A0 ^= pa ? delta : zero_block;  // set A0 as the label with signal bit 0
//     B0 ^= pb ? delta : zero_block;  // set B0 as the label with signal bit 0
//     bool barpa = pa ^ 1; //2 * barpa + \barpb is the position of 1 in the garbled table, in fact, bara and barb is the "a" and "b" in the paper "three halves makes a whole"
//     bool barpb = pb ^ 1; 

//     halfblock LA, RA, LB, RB, LC, RC, LDelta, RDelta;
//     splitBlock(A0, LA, RA);
//     splitBlock(B0, LB, RB);
//     splitBlock(delta, LDelta, RDelta); // split the input labels and the public offset

//     // PRG prg;
//     bool rand[2]; 
//     prg->random_bool(rand, 2); //choose the coefficient of generation R_$
//     rand[0] = rand[1] = 0;
//     bool VinvRp[30] = {
//         0, 0, 1, 0, 0, 0,
//         0, 1, 0, 0, 0, 0,
//         0, 0, 1, 0, 0, 1,
//         0, 1, 0, 0, 1, 0,
//         0, 0, 0, 0, 0, 0
//     };

//     bool VinvRa[30] = {
//         0, 0, 0, 0, 0, 0,
//         0, 0, 0, 0, 0, 0,
//         1, 1, 1, 0, 1, 1,
//         1, 0, 0, 1, 0, 1, 
//         0, 1, 1, 1, 1, 1
//     };

//     bool VinvRb[30] = {
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 0, 0, 
//         1, 0, 0, 1, 1, 0,
//         0, 1, 1, 1, 1, 1, 
//         1, 1, 1, 0, 1, 0 
//     };

//     bool VinvRr1[30] = {
//         1, 1, 1, 0, 0, 0, 
//         1, 0, 0, 1, 0, 0, 
//         0, 0, 0, 0, 0, 1, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 1, 0
//     };

//     bool VinvRr2[30] = {
//         1, 0, 0, 1, 0, 0, 
//         0, 1, 1, 1, 0, 0, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 1, 0, 
//         0, 0, 0, 0, 0, 1
//     };

//     bool VinvT[120] = {
//         0, 0, 0, 0, 1, 0,// VinvT00
//         0, 0, 0, 0, 0, 1, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 0, 0,

//         0, 0, 0, 0, 0, 0, // VinvT01 
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 0, 0,

//         0, 0, 0, 0, 0, 0, // VinvT10 
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 1, 1, 
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 1, 0,

//         0, 0, 0, 0, 0, 0, //VinvT11
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 0, 0, 
//         0, 0, 0, 0, 0, 0,
//         0, 0, 0, 0, 1, 0
//     };

//     // bool Vinv[40] = { // 5*8 matrix
//     // 1, 0, 0, 0, 0, 0, 0, 0,
//     // 0, 1, 0, 0, 0, 0, 0, 0,
//     // 1, 1, 0, 0, 1, 1, 0, 0,
//     // 1, 1, 1, 1, 0, 0, 0, 0,
//     // 0, 0, 0, 0, 1, 0, 1, 0};



//     // bool VinvM[30] = { // 5*6 matrix
//     // 1, 0, 0, 0, 1, 0,
//     // 0, 0, 1, 0, 1, 0,
//     // 1, 1, 0, 0, 0, 0,
//     // 0, 0, 1, 1, 0, 0,
//     // 0, 0, 0, 0, 1, 1};

    
//     bool VinvR[30] = {}; // 5*6 matrix

//     int pos = 30 * (2 * barpa + barpb);
//     for(int i = 0; i < 30; i ++)
//     {
//         VinvR[i] = (barpa & VinvRa[i]) ^ (barpb & VinvRb[i]) ^ (rand[0] & VinvRr1[i]) ^ (rand[1] & VinvRr2[i]) ^ VinvRp[i] ^ VinvT[pos + i];
//         std::cout << VinvR[i];
//     }
//     std::cout << std::endl;
//     // abort();
//     halfblock res1[5] = {};
//     halfblock res2[5] = {};
//     bool boolres1[5] = {};
//     bool boolres2[5] = {};
//     halfblock InputLabels[6] = {LA, RA, LB, RB, LDelta, RDelta};
//     gf2MatMul_bool_halfblock(res1, VinvR, InputLabels, 5, 6, 1); // compute Vinv \times R \times [A_0, B_0, \Delta]

//     block H[6];
// 	H[0] = A0;
// 	H[1] = A0 ^ delta;
// 	H[2] = B0;
// 	H[3] = B0 ^ delta;
//     H[4] = A0 ^ B0;
//     H[5] = A0 ^ B0 ^ delta;

// 	mitccrh->hash<3,2>(H); //compute H(A_0), H(A_1), H(B_0), H(B_1), H(A_0 \xor B_0), H(A_0 \xor B_1)
//     halfblock tmp;
//     halfblock HashLabels[6];
//     bool HF[6] = {};
//     for(int i = 0; i < 6; i++)
//     {
//         HF[i] = getMSB(H[i]);
//         HashLabels[i] = LBlock(H[i]); 
//         std::cout << "HashLabels" << HashLabels[i] << std::endl;
//     }
//     // splitBlock(H[0], HashLabels[0], tmp);
//     // splitBlock(H[1], HashLabels[1], tmp);
//     // splitBlock(H[2], HashLabels[2], tmp);
//     // splitBlock(H[3], HashLabels[3], tmp);
//     // splitBlock(H[4], HashLabels[4], tmp);
//     // splitBlock(H[5], HashLabels[5], tmp);


//     // gf2MatMul_bool_halfblock(res2, VinvM, HashLabels, 5, 6, 1);
//     res2[0] = HashLabels[0] ^ HashLabels[4];
//     res2[1] = HashLabels[2] ^ HashLabels[4];
//     res2[2] = HashLabels[0] ^ HashLabels[1];
//     res2[3] = HashLabels[2] ^ HashLabels[3];
//     res2[4] = HashLabels[4] ^ HashLabels[5];

//     boolres2[0] = HF[0] ^ HF[4];
//     boolres2[1] = HF[2] ^ HF[4];
//     boolres2[2] = HF[0] ^ HF[1];
//     boolres2[3] = HF[2] ^ HF[3];
//     boolres2[4] = HF[4] ^ HF[5];

//     // gf2MatMul_bool(boolres2, VinvM, HF, 5, 6, 1);
//     LC = res1[0] ^ res2[0];
//     RC = res1[1] ^ res2[1];
//     table[0] = res1[2] ^ res2[2];
//     table[1] = res1[3] ^ res2[3];
//     table[2] = res1[4] ^ res2[4];

//     // bool r[8] = {}; // the dicing bit
//     // r[0] = rand[0];
//     // r[1] = rand[1];
//     // r[2] = barpa ^ barpb ^ rand[0];
//     // r[3] = barpa ^ rand[1];
//     // r[4] = barpb ^ rand[0];
//     // r[5] = barpa ^ barpb ^ rand[1];
//     // r[6] = barpa ^ rand[0];
//     // r[7] = barpb ^ rand[1];
//     // gf2MatMul_bool(boolres2, Vinv, r, 5, 8, 1);

//     boolres1[0] = rand[0];
//     boolres1[1] = rand[1];
//     boolres1[2] = barpa;
//     boolres1[3] = barpb;
//     boolres1[4] = barpa^barpb;

//     for(int i = 0; i < 5; i++)
//     {
//         z[i] = boolres1[i] ^ boolres2[i];
//     }

//     block C = makeBlock(LC, RC);
//     std::cout << C << std::endl;
//     std::cout <<"LAbels" << std::endl << LA << std::endl << RA << std::endl << LB << std::endl << RB << std::endl << LDelta << std::endl << RDelta << std::endl;
//     abort();
// 	return C;
// }















template<typename T>
class TreGateGen:public CircuitExecution {
public:
	block delta;
	T * io;
	block constant[2];
    PRG prg;
	MITCCRH<24> mitccrh;
	TreGateGen(T * io) :io(io) {
		block tmp[2];
		PRG().random_block(tmp, 2);
		set_delta(tmp[0]);
		io->send_block(tmp+1, 1);
		mitccrh.setS(tmp[1]);
	}
	void set_delta(const block & _delta) {
		delta = set_bit(_delta, 0);
		PRG().random_block(constant, 2);
		io->send_block(constant, 2);
		constant[1] = constant[1] ^ delta;
	}
	block public_label(bool b) override {
		return constant[b];
	}
	block and_gate(const block& a, const block& b) override {
		halfblock table[3];
        bool z[5];
		block res = tregates_garble(a, b, delta, table, z, &mitccrh, &prg);
		io->send_halfblock(table, 3);
        io->send_bool(z, 5);
		return res;
	}
	block xor_gate(const block&a, const block& b) override {
		return a ^ b;
	}
	block not_gate(const block&a) override {
		return xor_gate(a, public_label(true));
	}
	uint64_t num_and() override {
		return mitccrh.gid/6;
	}
};
}
#endif// HALFGATE_GEN_H__





























































// inline block tregates_garble(block A0, block B0, block delta, halfblock *table, bool *z, MITCCRH<24> *mitccrh, PRG *prg) {
//     bool pa = getLSB(A0);
// 	bool pb = getLSB(B0);

//     A0 ^= pa ? delta : zero_block;  // set A0 as the label with signal bit 0
//     B0 ^= pb ? delta : zero_block;  // set B0 as the label with signal bit 0
//     bool barpa = pa ^ 1; //2 * barpa + \barpb is the position of 1 in the garbled table, in fact, bara and barb is the "a" and "b" in the paper "three halves makes a whole"
//     bool barpb = pb ^ 1; 

//     halfblock LA, RA, LB, RB, LC, RC, LDelta, RDelta;
//     splitBlock(A0, LA, RA);
//     splitBlock(B0, LB, RB);
//     splitBlock(delta, LDelta, RDelta); // split the input labels and the public offset


//     block H[6];
//     halfblock HashLabels[6];
//     halfblock tmp;
// 	H[0] = A0;
// 	H[1] = A0 ^ delta;
// 	H[2] = B0;
// 	H[3] = B0 ^ delta;
//     H[4] = B0 ^ 1;
//     H[5] = A0 ^ 1;

//     mitccrh->hash<3,2>(H);

//     splitBlock(H[0], HashLabels[0], tmp);
//     splitBlock(H[1], HashLabels[1], tmp);
//     splitBlock(H[2], HashLabels[2], tmp);
//     splitBlock(H[3], HashLabels[3], tmp);
//     splitBlock(H[4], HashLabels[4], tmp);
//     splitBlock(H[5], HashLabels[5], tmp);


//     table[0] = HashLabels[2] ^ HashLabels[3] ^ pa;
//     table[1] = HashLabels[0] ^ HashLabels[1] ;
//     table[2] = HashLabels[4] ^ HashLabels[5];

// 	return A0;
// }















