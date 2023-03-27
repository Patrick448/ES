#ifndef _LSODA_H_
#define _LSODA_H_
#include <stdlib.h>
#pragma once

#ifdef LIBLSODA_EXPORTS
#define LIBLSODA __declspec(dllexport)
#else
#define LIBLSODA __declspec(dllimport)
#endif

/* ************************************
 * 
 */
struct lsoda_opt_t {
	int ixpr;
	int mxstep;
	int mxhnil;
	int mxordn;
	int mxords;
	double tcrit;
	double h0;
	double hmax;
	double hmin;
	double hmxi;
	int itask;
	double *rtol;
	double *atol;
};

typedef int (*_lsoda_f) (double, double *, double *, double *);


struct lsoda_context_t {
	_lsoda_f function;
	double * data;
	int neq;
	int state;
	char * error;
/* private for lsoda */
	struct lsoda_common_t * common;
	struct lsoda_opt_t * opt;
};

#ifdef __cplusplus
extern "C" {
#endif

LIBLSODA int lsoda_prepare(struct lsoda_context_t * ctx, struct lsoda_opt_t * opt);
void lsoda_reset(struct lsoda_context_t * ctx);
LIBLSODA int lsoda(struct lsoda_context_t * ctx, double *y, double *t, double tout);
LIBLSODA void lsoda_free(struct lsoda_context_t * ctx);
void lsoda_free_opt(struct lsoda_opt_t * opt);

struct lsoda_context_t * lsoda_create_ctx();
struct lsoda_opt_t * lsoda_create_opt();

#ifdef __cplusplus
}
#endif

#endif
