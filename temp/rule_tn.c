#include <stddef.h>

typedef func;

#define MIN(a, b) ((a > b) ? (a) : (b))

double dot(double* a, double* b, size_t n);
double triple_dot(double* a, double* b, double* c, size_t n);

void matmul(double* out, double* a, double* b, size_t m, size_t n, size_t p);
void transpose(double* out, double* a, size_t m, size_t n);

void multiply(double* out, double* a, double* b, size_t n);

size_t sum(size_t* a, size_t n);

void rule_tn(
    double tune, size_t ndim, const double* vertex, double volume, size_t fdim, func integrand, int inkey, double* basval, double* rgnerr, size_t fvalt)
{
    size_t max_weights = 21;
    size_t max_rules = 7;
    size_t max_gen = 4;


    double* generators = (double*) calloc((max_gen + 1)*max_weights, sizeof(double));
    size_t* pts = (size_t*) calloc(max_weights, sizeof(size_t));

    size_t num_weights = 0;
    size_t oldn = 0;
    int oldkey = -1;
    int key;
    if (oldkey != inkey || oldn != ndim)
    {
        oldn = ndim;
        oldkey = inkey;
        if (inkey > 0 && inkey < 5)
            key = inkey;
        else    
            key = 3;
        
        double* weights = (double*) malloc(max_weights*max_rules*sizeof(double));
        size_t num_rules;
        ruleparams_tn(ndim, key, weights, generators, &num_weights, &num_rules, pts);

        double normcf = triple_dot(pts, weights, weights, num_weights);
        for (size_t k = 1; k < num_rules; ++k)
        {
            double* wt = (double*) malloc(max_weights*max_rules*sizeof(double));
            transpose(wt, weights, k - 2, max_weights);

            double* pts_w = (double*) malloc(max_weights*sizeof(double));
            for (size_t i = 0; i < max_weights; ++i)
                pts_w[i] = pts[i]*weights[i];

            double* alpha = (double*) malloc(max_rules*sizeof(double));
            matmul(&alpha[1], wt, pts_w, k - 1, max_weights, 1);

            double* temp = (double*) malloc(max_weights*sizeof(double));

            /* matmul here */

            for (size_t i = 0; i < max_weights; ++i)
                weights[k*max_weights + i] += temp[i];


            double normnl = triple_dot(pts, &weights[k*max_weights], &weights[k*max_weights], max_weights);
            for (size_t i = 0; i < max_weights; ++i)
                weights[k*max_weights + i] *= sqrt(normcf/normnl);
        }
        fvalt = sum(pts, num_weights);
    }

    double* rule = (double*) calloc(fdim*max_rules, sizeof(double));
    for (size_t k = 0; k < num_weights; ++k)
    {
        if (pts[k] > 0)
        {
            double* gtemp = (double*) malloc((ndim + 1)*sizeof(double));
            memcpy(gtemp, &generators[k*(max_gen + 1)], (MIN(ndim, max_gen - 1) + 1)*sizeof(double));
            if (ndim >= max_gen)
                for (size_t i = 0; i < ndim; ++i)
                    gtemp[max_gen + i] = generators[k*(max_gen + 1) + max_gen];
            symsmp_sum(basval, ndim, vertex, fdim, integrand, gtemp);
            for (size_t i = 0; i < fdim*max_rules; ++i)
                rule[i] += /* stuff here */0;
        }
    }
    memcpy(basval, rule, fdim*sizeof(double));
}