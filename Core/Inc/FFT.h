#ifndef INC_FFT_H_
#define INC_FFT_H_

#include "main.h"

typedef int8_t q8_t;
typedef int16_t q15_t;
typedef int32_t q31_t;

q15_t floatToQ15(float val);
float Q15ToFloat(q15_t val);

#define ADDQ15(x, y) __SSAT((q31_t)x + y, 16)
#define SUBQ15(x, y) __SSAT((q31_t)x - y, 16)
#define MULQ15(x, y) __SSAT(((q31_t)x * y) >> 15, 16)

#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...
#define  FT_DIRECT        -1    // Direct transform.
#define  FT_INVERSE        1    // Inverse transform.

void FFT(float* Rdat, float* Idat, uint16_t N, int8_t Ft_Flag);

uint16_t FFT_Q15(q15_t* Rdat, q15_t* Idat, uint16_t N, int8_t Ft_Flag);

#endif /* INC_FFT_H_ */
