/* problem3.c
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3, as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "problem3.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct problem {
    /*
     * IMPLEMENT HERE
     */
    double *matrix;
    int n; // number of nodes/groups
    int e; // number of components
};

struct solution {
    struct problem *prob;
    /*
     * IMPLEMENT HERE
     */
    int *nodes;
    int **groups;
    int *group_sizes;
    int cur_num_components;
    double objvalue;
    int evalv;    /* Flag indicating if the solution is evaluated */
    double objv;  /* Objective value */
    int evalLB;   /* Flag indicating if the lower bound is calculated */
    double objLB; /* Lower bound */
};

struct move {
    struct problem *prob;
    /*
     * IMPLEMENT HERE
     */
    int node;
    int group;
    int evalLBi[2]; /* Flag indicating if lower bound increment is evaluated for subneighbourhoods: { 0 - Add, 1 - Remove } */
    double objLBi;  /* Lower bound increment */
};

extern gsl_rng *rng; /* The single rng instance used by the whole code */

/*********************************/
/* ----- Utility functions ----- */
/*********************************/

/*
 * Return random integer in the range 0..N
 */
static int randint(const int n_max)
{
    return gsl_rng_uniform_int(rng, n_max + 1);
}

/*
 * Exchange the values of the ith and jth elements of an array
 */
static void swap(int *data, const int i, const int j)
{
    if (i == j)
        return;
    int val = data[i];
    data[i] = data[j];
    data[j] = val;
}

/*
 * Exchange the values of the ith and jth elements of a permutation.
 * Update the inverse permutation.
 */
static void swap_i(int *data, int *idata, const int i, const int j)
{
    if (i == j)
        return;
    swap(data, i, j);
    swap(idata, *(data + i), *(data + j));
}

int index_calc(int i, int j, int n) {
    if (i > j) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    int index = 0;
    for(int k = 1; k < i+1; ++k) {
        index += n - k;
    }
    index += j - i - 1;
    //return i*(n-1-i) + j-i;
    return index;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Problem instantiation
 */
struct problem *newProblem(const char *filename)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    int i = 0;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open file %s.\n", filename);
        return NULL;
    }
    int n = 0;
    if((getline(&line, &len, fp)) != -1) {
        n = atoi(line);
    }
    if (n == 0) {
        fprintf(stderr, "Invalid number of vertices: %d\n", n);
    }
    struct problem *p = NULL;
    p = (struct problem *) malloc(sizeof (struct problem));
    p->n = n;
    int matrix_size = n*(n-1)/2;
    p->matrix = (double *)malloc(matrix_size * sizeof(double));

    for (i = 1; i <= n; ++i) {
        p->e *= i;
    }

    i = 0;

    while ((getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
        char *token = strtok(line," ");
        int j = i;
        int value;
        while(token!=NULL)
        {
            if(j > i && j < n) {
                value = atoi(token);
                int index = index_calc(i, j, n);
                p->matrix[index] = (double)value;
            }
            token=strtok(NULL," ");
            ++j;
        }
        ++i;
    }
    fclose(fp);
    return p;
}

/**********************/
/* Problem inspection */
/**********************/

/*
 * Return the largest possible number of neighbours in a given subneighbourhood
 */
long getMaxNeighbourhoodSize(const struct problem *p, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
        return p->n;
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
        return 1;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getMaxNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*
 * Return the size of the ground set of the problem instance
 */
long getNumComponents(const struct problem *p)
{
    /*
     * IMPLEMENT HERE
     */
    return p->e;
}

/*
 * Return the largest number of components that a solution can potentially have.
 */
long getMaxSolutionSize(const struct problem *p)
{
    /*
     * IMPLEMENT HERE
     */
    return p->n;
}

/*********************/
/* Memory management */
/*********************/

/*
 * Allocate memory for a solution
 */
struct solution *allocSolution(struct problem *p)
{
    struct solution *s = malloc(sizeof(struct solution));
    s->prob = p;
    /*
     * IMPLEMENT HERE
     */
    return s;
}

/*
 * Allocate memory for a move
 */
struct move *allocMove(struct problem *p)
{
    struct move *v = malloc(sizeof(struct move));
    v->prob = p;
    /*
     * IMPLEMENT HERE
     */
    return v;
}

/*
 * Free the memory used by a problem
 */
void freeProblem(struct problem *p)
{
    /*
     * IMPLEMENT HERE
     */
    free(p);
}

/*
 * Free the memory used by a solution
 */
void freeSolution(struct solution *s)
{
    /*
     * IMPLEMENT HERE
     */
    free(s);
}

/*
 * Free the memory used by a move
 */
void freeMove(struct move *v)
{
    /*
     * IMPLEMENT HERE
     */
    free(v);
}

/*************/
/* Reporting */
/*************/

/*
 * Print user-formatted representation of problem instance
 */
void printProblem(const struct problem *p)
{
    /*
     * IMPLEMENT HERE
     */
}

/*
 * Print user-formatted representation of solution
 */
void printSolution(const struct solution *s)
{
    /*
     * IMPLEMENT HERE
     */
}

/*
 * Print user-formatted representation of move
 */
void printMove(const struct move *v)
{
    /*
     * IMPLEMENT HERE
     */
}

/***************************/
/* Operations on Solutions */
/***************************/

/*
 * Initialise empty solution
 */
struct solution *emptySolution(struct solution *s)
{
    /* solution s must have been allocated with allocSolution() */
    /*
     * IMPLEMENT HERE
     */
    return s;
}

/*
 * Copy the contents of the second argument to the first argument
 */
struct solution *copySolution(struct solution *dest, const struct solution *src)
{
    dest->prob = src->prob;
    /*
     * IMPLEMENT HERE
     */
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
    return dest;
}

/*
 * Solution evaluation
 */
double *getObjectiveVector(double *objv, struct solution *s)
{
    /* solution is unfeasible, cannot evaluate it */
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
        /*
         * IMPLEMENT HERE
         */
        *objv = s->objv = obj;
        s->evalv = 1;
    }
    return objv;
}

/*
 * Lower bound evaluation
 */
double *getObjectiveLB(double *objLB, struct solution *s)
{
    double obj = 0.0;
    if (s->evalLB) /* solution s is evaluated */
        *objLB = s->objLB;
    else { /* solution s is not evaluated */
        /*
         * IMPLEMENT HERE
         */
        *objLB = s->objLB = obj;
        s->evalLB = 1;
    }
    return objLB;
}

/*
 * Modify a solution in place by applying a move to it
 */
struct solution *applyMove(struct solution *s, const struct move *v, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
        i = 0;
        break;
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
        i = 1;
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to applyMove().\n");
        return NULL;
    }
    /* update state of evaluation */
    s->evalv = 0;
    if (s->evalLB && v->evalLBi[i])
        s->objLB += v->objLBi;
    else
        s->evalLB = 0;
    return s;
}

/*
 * Return true if a given solution is feasible or false if it is unfeasible
 */
int isFeasible(struct solution *s)
{
    /*
     * IMPLEMENT HERE
     */
}

/*
 * Reset the enumeration of a given subneighbourhood of a solution
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid neighbourhood passed to resetEnumMove().\n");
        break;
    }
    return NULL;
}

/*
 * Return the number of neighbours in a given subneighbouhood of a solution
 */
long getNeighbourhoodSize(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getNeighbourhoodSize().\n");
        break;
    }
    return -1;
}

/*
 * Enumerate the components of a solution that are in a given state
 */
long enumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid state passed to enumSolutionComponents().\n");
        break;
    }
    return -1;
}

/*
 * Reset the enumeration of the components of a solution that are in a given
 * state
 */
struct solution *resetEnumSolutionComponents(struct solution *s, const enum ComponentState st)
{
    switch (st) {
    case PRESENT:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid state passed to resetEnumSolutionComponents().\n");
        break;
    }
    return NULL;
}

/*
 * Heuristically constructs a feasible solution
 */
struct solution *heuristicSolution(struct solution *s)
{
    /*
     * IMPLEMENT HERE
     */
    return s;
}

/***********************/
/* Operations on Moves */
/***********************/

/*
 * Enumeration of a given subneighbourhood of a solution
 */
struct move *enumMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    /* subneighbourhood nh of solution is an empty set, cannot generate move */
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid neighbourhood passed to applyMove().\n");
        return NULL;
    }
    memset(v->evalLBi, 0, sizeof(int) * 2);
    return v;
}

/*
 * Copy the contents of the second argument to the first argument
 */
struct move *copyMove(struct move *dest, const struct move *src)
{
    dest->prob = src->prob;
    /*
     * IMPLEMENT HERE
     */
    memcpy(dest->evalLBi, src->evalLBi, 2 * sizeof(int));
    dest->objLBi = src->objLBi;
    return dest;
}

/*
 * Move evaluation
 */
double *getObjectiveLBIncrement(double *obji, struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    int i;
    switch (nh) {
    case ADD:
        i = 0;
        break;
    case REMOVE:
        i = 1;
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to getObjectiveLBIncrement().\n");
        return NULL;
    }
    if (v->evalLBi[i]) /* move v is evaluated */
        *obji = v->objLBi;
    else { /* move v is not evaluated */
        memset(v->evalLBi, 0, sizeof(int) * 2);
        switch (nh) {
        case ADD:
            /*
             * IMPLEMENT HERE
             */
            break;
        case REMOVE:
            /*
             * IMPLEMENT HERE
             */
            break;
        default:
            *obji = v->objLBi = DBL_MAX;
            break;
        }
        v->evalLBi[i] = 1;
    }
    return obji;
}

/*
 * Return the unique component identifier with respect to a given move
 */
long getComponentFromMove(const struct move *v)
{
    /*
     * IMPLEMENT HERE
     */
}

/*
 * Uniform random sampling with replacement of a given subneighbourhood of a
 * solution
 */
struct move *randomMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    /* subneighbourhood nh of solution is an empty set, cannot generate move */
    switch (nh) {
    case ADD:
        /*
         * IMPLEMENT HERE
         */
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
    default:
        fprintf(stderr, "Invalid neighbourhood passed to applyMove().\n");
        return NULL;
    }
    memset(v->evalLBi, 0, sizeof(int) * 2);
    return v;
}
