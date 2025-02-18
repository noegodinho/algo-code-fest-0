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
    int cur_num_groups;
    int cur_enumMoves;
    int cur_enumSolutionComponents;
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
    struct problem *p = (struct problem *) malloc(sizeof (struct problem));
    p->n = n;
    int matrix_size = n*(n-1)/2;
    p->matrix = (double *)malloc(matrix_size * sizeof(double));

    p->e = 0;

    int i;
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
    struct solution *s = (struct solution *)malloc(sizeof(struct solution));
    /*
     * IMPLEMENT HERE
     */
    int n;

    s->prob = p;
    n = p->n;

    s->nodes = (int *)malloc(n * sizeof(int));
    s->groups = (int **)malloc(n * sizeof(int *));
    s->group_sizes = (int *)malloc(n * sizeof(int));

    for(int i = 0; i < n; ++i) {
        s->groups[i] = (int *)malloc(n * sizeof(int));
    }

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
    free(p->matrix);
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
    int n = s->prob->n;
    free(s->nodes);

    for(int i = 0; i < n; ++i) {
        free(s->groups[i]);
    }
    
    free(s->groups);
    free(s->group_sizes);
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
    int n = p->n;
    printf("  \t");
    for(int i = 0; i < n; ++i) {
        printf( "%i\t", i);
    }
    printf("\n");
    for(int i = 0; i < n; ++i) {
        printf("%i|\t", i);
        for(int j = 0; j < n; ++j) {
            if (j > i) {
                printf("%i\t", (int)(p->matrix[index_calc(i, j, n)]));
            }
            else {
                printf(" \t");
            }
        }
        printf("\n");
    }
}

/*
 * Print user-formatted representation of solution
 */
void printSolution(const struct solution *s)
{
    /*
     * IMPLEMENT HERE
     */
    int n = s->prob->n;
    printf("Node: \t");
    for(int i = 0; i < n; ++i) {
        printf("%i\t", i);
    }
    printf("\n");
    printf("Group:\t");
    for(int i = 0; i < n; ++i) {
        printf("%i\t", s->nodes[i]);
    }
    printf("\n\n");
    for(int i = 0; i < n; ++i) {
        printf("Group %i: ", i);
        for(int j = 0; j < s->group_sizes[i]; ++j) {
            printf("%i ", s->groups[i][j]);
        }
        printf("\n");
    }
}

/*
 * Print user-formatted representation of move
 */
void printMove(const struct move *v)
{
    /*
     * IMPLEMENT HERE
     */
    printf("Move (node, group): %d, %d\n", v->node, v->group);
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
    for(int i = 0; i < s->prob->n; ++i){
        s->nodes[i] = -1;
        s->group_sizes[i] = 0;
    }

    s->cur_num_components = 0;
    s->cur_num_groups = 0;
    s->objvalue = 0;
    s->evalv = 0;
    s->objv = 0;
    s->evalLB = 0;
    s->objLB = 0;
    s->cur_enumMoves = 0;
    s->cur_enumSolutionComponents = 0;

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
    memcpy(dest->nodes, src->nodes, src->prob->n * sizeof (int));

    for(int i = 0; i < src->prob->n; ++i) {
        memcpy(dest->groups[i], src->groups[i], src->group_sizes[i] * sizeof (int));
        dest->group_sizes[i] = src->group_sizes[i];
    }

    dest->cur_num_components = src->cur_num_components;
    dest->cur_num_groups = src->cur_num_groups;
    dest->objvalue = src->objvalue;
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
    dest->cur_enumMoves = src->cur_enumMoves;
    dest->cur_enumSolutionComponents = src->cur_enumSolutionComponents;
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
        for(int i = 0; i < s->prob->n; ++i) {
            double groupObjVal = 0;
            for(int j = 0; j < s->group_sizes[i]; ++j) {
                for(int k = j+1; k < s->group_sizes[i]; ++k) {
                    int index = index_calc(s->groups[i][j], s->groups[i][k], s->prob->n);
                    groupObjVal += s->prob->matrix[index];
                }
            }
            obj += groupObjVal;
        }
        *objv = s->objv = -obj;
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
        double value;
        for(int j = s->cur_num_components; j < s->prob->n; ++j) {
            for(int k = 0; k < s->prob->n; ++k) {
                int index = index_calc(j, k, s->prob->n);
                value = s->prob->matrix[index];
                if (value > 0) {
                    obj += value;
                }
            }
        }
        *objLB = s->objLB = -obj;
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
        s->nodes[v->node] = v->group;       
        if(s->group_sizes[v->group] > s->prob->n){
            printf("This should not happen.\n");
        }
        s->groups[v->group][s->group_sizes[v->group]++] = v->node;
        i = 0;
        break;
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
        --s->group_sizes[v->group];
        s->nodes[v->node] = -1;
        for(int x = 0; x <= s->group_sizes[v->group]; x++){
            if(s->groups[v->group][x] == v->node){
                s->groups[v->group][x] = s->groups[v->group][s->group_sizes[v->group]];
                break;
            }
        }
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

    s->cur_enumMoves = 0;
    s->cur_enumSolutionComponents = 0;
    return s;
}

/*
 * Return true if a given solution is feasible or false if it is unfeasible
 */
int isFeasible(struct solution *s)
{
    return s->cur_num_components == s->prob->n;
}

/*
 * Reset the enumeration of a given subneighbourhood of a solution
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        //s->cur_enumMoves = 0;
        //return s;
    case REMOVE:
        /*
         * IMPLEMENT HERE
         */
        s->cur_enumMoves = 0;
        return s;
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
        if (s->cur_num_components == s->prob->n) {
            return 0;
        }
        else {
            return s->cur_num_groups + 1;
        }
    case REMOVE:
        if (s->cur_num_components == 0) {
            return 0;
        }
        else {
            return 1;
        }
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
        if (s->cur_num_components == 0) {
            return -1;
        }
        else {
            int identifier = index_calc(s->nodes[s->cur_enumSolutionComponents],s->cur_enumSolutionComponents, s->prob->n);
            s->cur_enumSolutionComponents++;
            return identifier;
        } 
    default:
        fprintf(stderr, "Invalid state passed to enumSolutionComponents().\n");
        break;
    }
    return -1;
}

long getComponentFromMove(const struct move *v) 
{
    return index_calc(v->group, v->node, v->prob->n);
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
        s->cur_enumSolutionComponents = 0;
        return s;
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
    struct move  *v = allocMove(s->prob);
    if (s->cur_num_components < s->prob->n) {
        for(int i = s->cur_num_components; i < s->prob->n ; i++) {
            v->node = i;
            v->group = randint(s->cur_num_groups);
            s->nodes[s->cur_num_components] = v->group;
            if (v->group == s->cur_num_components) {
                s->cur_num_groups++;
            }
            s->cur_num_components++;
            s->groups[v->group][s->group_sizes[v->group]++] = v->node;
            s->evalv = 0;
        }
    }

    s->cur_enumMoves = 0;
    s->cur_enumSolutionComponents = 0;
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
        if (s->cur_enumMoves <= s->cur_num_groups && s->cur_num_components < s->prob->n) {
            v->node = s->cur_num_components;
            v->group = s->cur_enumMoves;
            s->cur_enumMoves++;
        }
        else {
            return NULL;
        }

        break;
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
    dest->node = src->node;
    dest->group = src->group;
    dest->evalLBi[0] = src->evalLBi[0];
    dest->evalLBi[1] = src->evalLBi[1];
    dest->objLBi = src->objLBi;
    return dest;
}

/*
 * Move evaluation
 */
static double lbi_add(const struct move *v, const struct solution *s)
{
    int node = v->node;
    int group = v->group;
    int index;
    int val;
    double obj = 0.0;

    for(int i = 0; i < s->group_sizes[group]; ++i){
        if(s->groups[group][i] != node){
            index = index_calc(s->groups[group][i], node, s->prob->n);
            val = s->prob->matrix[index];

            if(val < 0){
                obj += val;
            }
        }
    }

    for(int i = 0; i < s->prob->n; ++i){
        if(i != node && s->nodes[i] != group){
            index = index_calc(i, node, s->prob->n);
            val = s->prob->matrix[index];

            if(val > 0){
                obj -= val;
            }
        }
    }

    return obj;
}

static double lbi_remove(const struct move *v, const struct solution *s)
{
    int node = v->node;
    int group = v->group;
    int index;
    int val;
    double obj = 0.0;

    for(int i = 0; i < s->group_sizes[group]; ++i){
        if(s->groups[group][i] != node){
            index = index_calc(s->groups[group][i], node, s->prob->n);
            val = s->prob->matrix[index];

            if(val > 0){
                obj += val;
            }
        }
        
    }

    for(int i = 0; i < s->prob->n; ++i){
        if(i != node && s->nodes[i] != group){
            index = index_calc(i, node, s->prob->n);
            val = s->prob->matrix[index];

            if(val > 0){
                obj -= val;
            }
        }
    }

    return obj;
}

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
            *obji = v->objLBi = lbi_add(v, s);
            break;
        case REMOVE:
            /*
             * IMPLEMENT HERE
             */
            *obji = v->objLBi = lbi_remove(v, s);
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
 * Uniform random sampling with replacement of a given subneighbourhood of a
 * solution
 */
struct move *randomMove(struct move *v, struct solution *s, const enum SubNeighbourhood nh)
{
    /* subneighbourhood nh of solution is an empty set, cannot generate move */
    switch (nh) {
    case ADD:
        if (s->cur_num_components < s->prob->n) {
            v->node = s->cur_num_components;
            v->group = randint(s->cur_num_groups);
        }
        else {
            return NULL;
        }
        break;
    case REMOVE:
        if (s->cur_num_components > 0) {
            v->node = s->cur_num_components - 1;
            v->group = s->nodes[s->cur_num_components - 1];
        }
        else {
            return NULL;
        }
        break;
    default:
        fprintf(stderr, "Invalid neighbourhood passed to applyMove().\n");
        return NULL;
    }
    memset(v->evalLBi, 0, sizeof(int) * 2);
    return v;
}
