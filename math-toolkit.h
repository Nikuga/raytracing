#ifndef __RAY_MATH_TOOLKIT_H
#define __RAY_MATH_TOOLKIT_H

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <immintrin.h>	

static inline __forceinline
void normalize(double *v)
{
    double d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    assert(d != 0.0 && "Error calculating normal");

    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}

static inline __forceinline
double length(const double *v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static inline __forceinline
void add_vector(const double *a, const double *b, double *out)
{
    //for (int i = 0; i < 3; i++)
    //    out[i] = a[i] + b[i];
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
}

static inline __forceinline
void subtract_vector(const double *a, const double *b, double *out)
{
    //for (int i = 0; i < 3; i++)
    //    out[i] = a[i] - b[i];
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
   
}

static inline __forceinline
void multiply_vectors(const double *a, const double *b, double *out)
{
    //for (int i = 0; i < 3; i++)
    //    out[i] = a[i] * b[i];
    out[0] = a[0] * b[0];
    out[1] = a[1] * b[1];
    out[2] = a[2] * b[2];

}

static inline __forceinline
void multiply_vector(const double *a, double b, double *out)
{
    //for (int i = 0; i < 3; i++)
    //    out[i] = a[i] * b;
    out[0] = a[0] * b;
    out[1] = a[1] * b;
    out[2] = a[2] * b;
}

static inline __forceinline
void cross_product(const double *v1, const double *v2, double *out)
{	
    out[0] = v1[1] * v2[2] - v1[2] * v2[1];
    out[1] = v1[2] * v2[0] - v1[0] * v2[2];
    out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static inline __forceinline
double dot_product(const double *v1, const double *v2)
{
    double dp = 0.0;
    //for (int i = 0; i < 3; i++)
    //    dp += v1[i] * v2[i];
    dp += v1[0] * v2[0];
    dp += v1[1] * v2[1];
    dp += v1[2] * v2[2];
	return dp;
	
	/*
	// method 1
	__m256d x = _mm256_set_pd( v1[0], v1[1], v1[2], 0.0);
	__m256d y = _mm256_set_pd( v2[0], v2[1], v2[2], 0.0);
    __m256d xy = _mm256_mul_pd(x, y);
    __m256d temp = _mm256_hadd_pd(xy, xy);
    __m128d hi128 = _mm256_extractf128_pd(temp, 1);
    __m128d low128 = _mm256_extractf128_pd(temp, 0);
    __m128d dotproduct = _mm_add_pd(low128, hi128);
	double dp = _mm_cvtsd_f64(dotproduct);
	return dp;
	*/
	
	/*
	//method 2
	double tmp[4] __attribute__((aligned(32)));
	__m256d x = _mm256_set_pd( v1[0], v1[1], v1[2], 0.0);
	__m256d y = _mm256_set_pd( v2[0], v2[1], v2[2], 0.0);
    __m256d xy = _mm256_mul_pd(x, y);
    _mm256_store_pd(tmp, xy);
    //dp = tmp[0] + tmp[1] + tmp[2] + tmp[3];
	
    return (tmp[0] + tmp[1] + tmp[2] + tmp[3]);
	*/
	
	/*
	//method 3
	double tmp[4] __attribute__((aligned(32)));
	__m256d x = _mm256_set_pd( v1[0], v1[1], v1[2], 0.0);
	__m256d y = _mm256_set_pd( v2[0], v2[1], v2[2], 0.0);
    __m256d xy = _mm256_mul_pd(x, y);	
	__m256d lvTemp1 = _mm256_permute_pd(xy,5);
	__m256d lvTemp2 = _mm256_permute2f128_pd(lvTemp1,lvTemp1,1);
	lvTemp2 = _mm256_permute_pd(lvTemp2,5);
	__m256d sum = _mm256_add_pd(xy,_mm256_add_pd(lvTemp1, lvTemp2)); 
    _mm256_store_pd(tmp, sum);
	
    return (tmp[3]);
	*/
	
}

static inline __forceinline
void scalar_triple_product(const double *u, const double *v, const double *w,
                           double *out)
{
    cross_product(v, w, out);
    multiply_vectors(u, out, out);
}

static inline __forceinline
double scalar_triple(const double *u, const double *v, const double *w)
{
    double tmp[3];
    cross_product(w, u, tmp);
    return dot_product(v, tmp);
}

#endif
