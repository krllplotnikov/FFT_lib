#include "FFT.h"

void FFT(float *Rdat, float *Idat, uint16_t N, int8_t Ft_Flag) {
	if ((Rdat == NULL) || (Idat == NULL))
		return;
	if ((N > 16384) || (N < 1))
		return;
	if (!NUMBER_IS_2_POW_K(N))
		return;
	if ((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE))
		return;

	uint8_t LogN;
	uint32_t i, j, n, k, io, ie, in, nn;
	float ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

	static const float Rcoef[14] = { -1.0000000000000000F, 0.0000000000000000F,
			0.7071067811865475F, 0.9238795325112867F, 0.9807852804032304F,
			0.9951847266721969F, 0.9987954562051724F, 0.9996988186962042F,
			0.9999247018391445F, 0.9999811752826011F, 0.9999952938095761F,
			0.9999988234517018F, 0.9999997058628822F, 0.9999999264657178F };
	static const float Icoef[14] = { 0.0000000000000000F, -1.0000000000000000F,
			-0.7071067811865474F, -0.3826834323650897F, -0.1950903220161282F,
			-0.0980171403295606F, -0.0490676743274180F, -0.0245412285229122F,
			-0.0122715382857199F, -0.0061358846491544F, -0.0030679567629659F,
			-0.0015339801862847F, -0.0007669903187427F, -0.0003834951875714F };

	uint8_t arrLogN[14] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
	uint16_t arrN[14] = { 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
			8192, 16384 };

	for (uint8_t i = 0; i < 14; i++) {
		if (N == arrN[i])
			LogN = arrLogN[i];
	}

	nn = N >> 1;
	ie = N;
	for (n = 1; n <= LogN; n++) {
		rw = Rcoef[LogN - n];
		iw = Icoef[LogN - n];
		if (Ft_Flag == FT_INVERSE)
			iw = -iw;
		in = ie >> 1;
		ru = 1.0F;
		iu = 0.0F;
		for (j = 0; j < in; j++) {
			for (i = j; i < N; i += ie) {
				io = i + in;
				rtp = Rdat[i] + Rdat[io];
				itp = Idat[i] + Idat[io];
				rtq = Rdat[i] - Rdat[io];
				itq = Idat[i] - Idat[io];
				Rdat[io] = rtq * ru - itq * iu;
				Idat[io] = itq * ru + rtq * iu;
				Rdat[i] = rtp;
				Idat[i] = itp;
			}

			sr = ru;
			ru = ru * rw - iu * iw;
			iu = iu * rw + sr * iw;
		}

		ie >>= 1;
	}

	for (j = i = 1; i < N; i++) {
		if (i < j) {
			io = i - 1;
			in = j - 1;
			rtp = Rdat[in];
			itp = Idat[in];
			Rdat[in] = Rdat[io];
			Idat[in] = Idat[io];
			Rdat[io] = rtp;
			Idat[io] = itp;
		}

		k = nn;

		while (k < j) {
			j = j - k;
			k >>= 1;
		}

		j = j + k;
	}

	if (Ft_Flag == FT_DIRECT)
		return;

	rw = 1.0F / N;

	for (i = 0; i < N; i++) {
		Rdat[i] *= rw;
		Idat[i] *= rw;
	}

	return;
}

uint16_t FFT_Q15(q15_t *Rdat, q15_t *Idat, uint16_t N, int8_t Ft_Flag) {
	if ((Rdat == NULL) || (Idat == NULL))
		return 0;
	if ((N > 16384) || (N < 1))
		return 0;
	if (!NUMBER_IS_2_POW_K(N))
		return 0;
	if ((Ft_Flag != FT_DIRECT) && (Ft_Flag != FT_INVERSE))
		return 0;

	uint8_t LogN;
	q15_t i, j, n, k, io, ie, in, nn;
	q15_t ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

	static const q15_t Rcoef[14] = { -32768, 0, 23170, 30273, 32138, 32610,
			32728, 32758, 32765, 32767, 32767, 32767, 32767, 32767 };
	static const q15_t Icoef[14] = { 0, -32768, -23170, -12539, -6392, -3211,
			-1607, -804, -402, -201, -100, -50, -25, -12 };

	uint8_t arrLogN[14] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };
	uint16_t arrN[14] = { 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
			8192, 16384 };

	for (uint8_t i = 0; i < 14; i++) {
		if (N == arrN[i])
			LogN = arrLogN[i];
	}

	nn = N >> 1;
	ie = N;
	for (n = 1; n <= LogN; n++) {
		if (Ft_Flag == FT_DIRECT){
			for(uint16_t i = 0; i < N; i++){
				Rdat[i] >>= 1;
				Idat[i] >>= 1;
			}
		}
		rw = Rcoef[LogN - n];
		iw = Icoef[LogN - n];
		if (Ft_Flag == FT_INVERSE)
			iw = -iw;
		in = ie >> 1;
		ru = 32767;
		iu = 0;
		for (j = 0; j < in; j++) {
			for (i = j; i < N; i += ie) {
				io = i + in;
				rtp = ADDQ15(Rdat[i], Rdat[io]);
				itp = ADDQ15(Idat[i], Idat[io]);
				rtq = SUBQ15(Rdat[i], Rdat[io]);
				itq = SUBQ15(Idat[i], Idat[io]);
				Rdat[io] = SUBQ15(MULQ15(rtq, ru), MULQ15(itq, iu));
				Idat[io] = ADDQ15(MULQ15(itq, ru), MULQ15(rtq, iu));
				Rdat[i] = rtp;
				Idat[i] = itp;
			}

			sr = ru;
			ru = SUBQ15(MULQ15(ru, rw), MULQ15(iu, iw));
			iu = ADDQ15(MULQ15(iu, rw), MULQ15(sr, iw));
		}

		ie >>= 1;
	}

	for (j = i = 1; i < N; i++) {
		if (i < j) {
			io = i - 1;
			in = j - 1;
			rtp = Rdat[in];
			itp = Idat[in];
			Rdat[in] = Rdat[io];
			Idat[in] = Idat[io];
			Rdat[io] = rtp;
			Idat[io] = itp;
		}

		k = nn;

		while (k < j) {
			j = j - k;
			k >>= 1;
		}

		j = j + k;
	}

	if (Ft_Flag == FT_DIRECT)
		return N;

	rw = 32767 / N;

	for (i = 0; i < N; i++) {
		Rdat[i] = MULQ15(Rdat[i], rw);
		Idat[i] = MULQ15(Idat[i], rw);
	}

	return 0;
}


q15_t floatToQ15(float val){
	if (val < -1.0f)
	{
		return (q15_t)0x8000;
	}
	else if (val >= 1.0f)
	{
		return (q15_t)0x7FFF;
	}
	else
	{
		return (q15_t)(val * 32768.0f); // val * 2^16
	}
}

float Q15ToFloat(q15_t val){
	return ((float)val * 3.0518043793e-05f);
}
