#ifndef EMP_MITCCRH_H__
#define EMP_MITCCRH_H__
#include "emp-tool/utils/aes_opt.h"
#include <stdio.h>

namespace emp {

/*
 * [REF] Implementation of "Better Concrete Security for Half-Gates Garbling (in the Multi-Instance Setting)"
 * https://eprint.iacr.org/2019/1168.pdf
 */

template<int BatchSize = 8>
class MITCCRH { public:
	AES_KEY scheduled_key[BatchSize];
	block keys[BatchSize];
	int key_used = BatchSize;
	block start_point;
	uint64_t gid = 0;

	void setS(block sin) {
		this->start_point = sin;
	}

	void renew_ks(uint64_t gid) {
		this->gid = gid;
		renew_ks();
	}

	void renew_ks() {
		for(int i = 0; i < BatchSize; ++i)
			keys[i] = start_point ^ makeBlock(gid++, 0);
		AES_opt_key_schedule<BatchSize>(keys, scheduled_key);
		key_used = 0;
	}

	template<int K, int H>
	void hash_cir(block * blks) {
		for(int i = 0; i < K*H; ++i)
			blks[i] = sigma(blks[i]);
		hash<K, H>(blks);
	}

	template<int K, int H>
	void hash(block * blks) {
		assert(K <= BatchSize);
		assert(BatchSize % K == 0);
		if(key_used == BatchSize) renew_ks();

		block tmp[K*H];
		for(int i = 0; i < K*H; ++i)
			tmp[i] = blks[i];
		
		ParaEnc<K,H>(tmp, scheduled_key+key_used);
		key_used += K;
		
		for(int i = 0; i < K*H; ++i)
			blks[i] = blks[i] ^ tmp[i];
	}








/**
	 * AesniTmmoBatch 统一模板函数 (固定密钥版)
	 * * @tparam NumTweaks       本次调用涉及的不同 Tweak 的数量 (例如 3)
	 * @tparam BlocksPerTweak  每个 Tweak 应用于多少个数据块 (Batch3则为1, Batch6则为2)
	 * * 推导属性: TotalBlocks (总处理Block数) = NumTweaks * BlocksPerTweak
	 */
	template <size_t NumTweaks, size_t BlocksPerTweak>
	void AesniTmmoBatch_FixedKey(void* input) {
		// 计算总处理的数据块数量 (编译期常量)
		constexpr size_t TotalBlocks = NumTweaks * BlocksPerTweak;
		gid += NumTweaks;
		// -------------------------------------------------------------------------
		// 1. 静态常量：固定的轮密钥 (Fixed Round Keys)
		// -------------------------------------------------------------------------
		alignas(16) static const std::array<__m128i, 11> round_keys = []() {
			std::array<__m128i, 11> keys;
			
			// [Step 1] 主密钥设为全 0
			keys[0] = _mm_setzero_si128(); 

			// [Step 2] 密钥扩展 Lambda
			auto expand = [](__m128i k, __m128i key_gen) {
				key_gen = _mm_shuffle_epi32(key_gen, _MM_SHUFFLE(3, 3, 3, 3));
				k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
				k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
				k = _mm_xor_si128(k, _mm_slli_si128(k, 4));
				return _mm_xor_si128(k, key_gen);
			};

			// [Step 3] 展开的密钥扩展
			keys[1]  = expand(keys[0], _mm_aeskeygenassist_si128(keys[0], 0x01));
			keys[2]  = expand(keys[1], _mm_aeskeygenassist_si128(keys[1], 0x02));
			keys[3]  = expand(keys[2], _mm_aeskeygenassist_si128(keys[2], 0x04));
			keys[4]  = expand(keys[3], _mm_aeskeygenassist_si128(keys[3], 0x08));
			keys[5]  = expand(keys[4], _mm_aeskeygenassist_si128(keys[4], 0x10));
			keys[6]  = expand(keys[5], _mm_aeskeygenassist_si128(keys[5], 0x20));
			keys[7]  = expand(keys[6], _mm_aeskeygenassist_si128(keys[6], 0x40));
			keys[8]  = expand(keys[7], _mm_aeskeygenassist_si128(keys[7], 0x80));
			keys[9]  = expand(keys[8], _mm_aeskeygenassist_si128(keys[8], 0x1b));
			keys[10] = expand(keys[9], _mm_aeskeygenassist_si128(keys[9], 0x36));

			return keys;
		}(); 

		// -------------------------------------------------------------------------
		// 2. 准备数据
		// -------------------------------------------------------------------------
		alignas(16) std::array<__m128i, TotalBlocks> wb_1;
		alignas(16) std::array<__m128i, TotalBlocks> wb_2;
		auto input_pointer = reinterpret_cast<__m128i*>(input);
		
		// 计算基础 Tweak
		// 这一批消耗了 NumTweaks 个索引
		uint64_t base_tweak = this->gid * NumTweaks;

		// -------------------------------------------------------------------------
		// 3. 第一层 AES: wb_1 <- \pi(x)
		// -------------------------------------------------------------------------
		// 编译器会自动展开这些循环 (Full Loop Unrolling)
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_xor_si128(input_pointer[j], round_keys[0]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[1]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[2]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[3]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[4]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[5]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[6]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[7]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[8]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenc_si128(wb_1[j], round_keys[9]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_1[j] = _mm_aesenclast_si128(wb_1[j], round_keys[10]);

		// -------------------------------------------------------------------------
		// 4. 应用 Tweak: wb_2 <- wb_1 ^ Tweak
		// -------------------------------------------------------------------------
		// 外层循环：遍历每一个 Tweak
		for (size_t t = 0; t < NumTweaks; ++t) {
			// 构建当前 Tweak (高位为0, 低位为 base + t)
			__m128i tweak_block = _mm_set_epi64x(0, base_tweak + t);
			
			// 内层循环：将当前 Tweak 应用到对应的几个 block 上
			for (size_t k = 0; k < BlocksPerTweak; ++k) {
				size_t idx = t * BlocksPerTweak + k; // 计算扁平化索引
				wb_2[idx] = _mm_xor_si128(wb_1[idx], tweak_block);
			}
		}

		// -------------------------------------------------------------------------
		// 5. 第二层 AES: wb_2 <- \pi(wb_2) 并输出
		// -------------------------------------------------------------------------
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_xor_si128(wb_2[j], round_keys[0]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[1]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[2]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[3]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[4]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[5]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[6]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[7]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[8]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenc_si128(wb_2[j], round_keys[9]);
		for (size_t j = 0; j < TotalBlocks; ++j) wb_2[j] = _mm_aesenclast_si128(wb_2[j], round_keys[10]);

		// store result
		for (size_t j = 0; j < TotalBlocks; ++j) input_pointer[j] = _mm_xor_si128(wb_2[j], wb_1[j]);
	}


























	
};
}
#endif// MITCCRH_H__
