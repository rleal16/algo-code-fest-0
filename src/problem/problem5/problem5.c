/* problem5.c
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

#include "problem5.h"
#include <float.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX(X, Y) ((X) > (Y)) ? (X) : (Y)

struct problem {
  int r, c, f, n, b, t;
  int **rides;
};

struct solution {
  struct problem *prob;
  int *rides;
  int *cars;
  int n_cars; // the number of cars present in the solution
  int n_rides;
  int last_eval_ride;
  int evalv;    /* Flag indicating if the solution is evaluated */
  double objv;  /* Objective value */
  int evalLB;   /* Flag indicating if the lower bound is calculated */
  double objLB; /* Lower bound */
  int addMoveCounter;     // used by enumMove()
  int removeMoveCounter;  // used by enumMove()
  int enumComponentState; /* number of components left to enumerate */
};

struct move {
  struct problem *prob;
  int componentId;
  int previousRide; /* 0..N-1 if previous is a ride or N if previous is a new car (newRide is the first ride in the car) */
  int newRide;
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

static void updateCarPositions(struct solution *s){
  struct problem *p = s->prob;
  int cars = 0;
  int n = p->f+p->n;
  for(int i = 0; i<n; i++){
    if(s->rides[i] >= p->n){
        s->cars[cars] = i;
        cars++;
    }
  }
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

static int pairing(int x, int y) {
    return x > y ? x * (x - 1) / 2 + y : y * (y - 1) / 2 + x;
}

/*************************/
/* Problem instantiation */
/*************************/

/*
 * Problem instantiation
 */
struct problem *newProblem(const char *filename)
{
    FILE *infile;
    int i;
    struct problem *p = NULL;
    infile = fopen(filename, "r");

    if (infile) {
        struct problem *p = (struct problem*) malloc(sizeof(struct problem));

        if(fscanf(infile, "%d %d %d %d %d %d", &(p->r), &(p->c), &(p->f), &(p->n), &(p->b), &(p->t)) != 6)
          return NULL;
        p->rides = (int**) malloc(sizeof(int*)*(p->n));

        //p->rides = malloc(sizeof(int)*(p->n));
        for(i = 0; i < p->n; i++){
          p->rides[i] = malloc(sizeof(int*)*6);
          if(fscanf(infile, "%d %d %d %d %d %d",
          &(p->rides[i][0]), &(p->rides[i][1]), &(p->rides[i][2]), &(p->rides[i][3]), &(p->rides[i][4]), &(p->rides[i][5])) != 6)
            return NULL;
        }

        fclose(infile);

        return p;
    }
    else
        fprintf(stderr, "Cannot open file %s.\n", filename);
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
        return p->n + p->n;
    case REMOVE:
        if(s->n_cars==1 && s->n_rides==0)
            return 0;
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
    return (p->n)*(1+p->n)/2;
}

/*
 * Return the largest number of components that a solution can potentially have.
 */
long getMaxSolutionSize(const struct problem *p)
{
    return (p->n);
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
    s->rides = malloc(sizeof(int)*(p->n + p->f));
    s->cars = malloc(sizeof(int)*(p->f));
    return s;
}

/*
 * Allocate memory for a move
 */
struct move *allocMove(struct problem *p)
{
    struct move *v = malloc(sizeof(struct move));
    v->prob = p;
    return v;
}

/*
 * Free the memory used by a problem
 */
void freeProblem(struct problem *p)
{
    int i;
    for(i = 0; i < p->n; i++){
        free(p->rides[i]);
    }
    free(p->rides);
    free(p);

}

/*
 * Free the memory used by a solution
 */
void freeSolution(struct solution *s)
{
  free(s->rides);
  free(s->cars);
  free(s);
}

/*
 * Free the memory used by a move
 */
void freeMove(struct move *v)
{
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
  int i;
  printf("%d %d %d %d %d %d\n", p->r, p->c, p->f, p->n, p->b, p->t);
  for(i = 0; i < p->n; i++){
    printf( "%d %d %d %d %d %d\n",
    p->rides[i][0], p->rides[i][1], p->rides[i][2], p->rides[i][3], p->rides[i][4], p->rides[i][5]);
  }
}

/*
 * Print user-formatted representation of solution
 */
void printSolution(const struct solution *s)
{
    for (int i = 0; i < s->n_cars + s->n_rides; ++i) {
        if (s->rides[i] >= s->prob->n)
            printf("\nthis vehicle has the following rides assigned:");
        else
            printf(" %d", s->rides[i]);
        if (s->evalv)
            printf("\nobjv = %.1lf", s->objv);
        if (s->evalLB)
            printf("\nobjLB = %.1lf", s->objv);
        printf("\n\n");
    }
}

/*
 * Print user-formatted representation of move
 */
void printMove(const struct move *v)
{
    if (v->previousRide < v->prob->n)
        printf("New ride %d after ride %d\n", v->newRide, v->previousRide);
    else
        printf("New car starting with ride %d\n", v->newRide);
    printf("component identifier %d\n", v->componentId);
    printf("increment %1.lf ", v->objLBi);
    if (v->evalLBi[0])
        printf("ADD\n");
    else
        printf("REMOVE\n");

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
  struct problem *p = s->prob;
  for(int i = 0; i < (p->n + p->f); i++){
    s->rides[i] = i;
  }
  s->rides[0] = p->n;
  s->rides[p->n] = 0;



  updateCarPositions(s);
  s->n_cars = 1;
  s->n_rides = 0;
  /** For Testing **/
  s->n_rides+=20;
  s->evalv = 0;
  s->evalLB = 0;
  s->addMoveCounter = 0;
  s->removeMoveCounter = 0;
  s->enumComponentState = s->n_rides;
  return s;
}

/*
 * Copy the contents of the second argument to the first argument
 */
struct solution *copySolution(struct solution *dest, const struct solution *src)
{
    dest->prob = src->prob;
    memcpy(dest->rides, src->rides, sizeof(int) * (src->prob->n + src->prob->f));
    memcpy(dest->cars, src->cars, sizeof(int) * (src->prob->f));
    dest->n_cars = src->n_cars;
    dest->n_rides = src->n_rides;
    dest->last_eval_ride = src->last_eval_ride;
    dest->evalv = src->evalv;
    dest->objv = src->objv;
    dest->evalLB = src->evalLB;
    dest->objLB = src->objLB;
    return dest;
}

/*
 * Solution evaluation
 */
 static double evaluateRide(int *step, int i, int *vh_pos, struct solution *s){
   struct problem *p = s->prob;
   double int_score = 0.0;
   int **rides = p->rides;
   int dist;
   int ride = s->rides[i];

   (*step) += abs(rides[ride][0]-vh_pos[0]) + abs(rides[ride][1]-vh_pos[1]); //Go to the start intersection

   //Give bonus if ride has started on time
   if((*step) <= rides[ride][4]){
     int_score+=p->b;
   }

   (*step) += MAX(0, rides[ride][4] - (*step)); //Wait for earliest start

   dist = abs(rides[ride][0]-rides[ride][2]) + abs(rides[ride][1]-rides[ride][3]);
   (*step) += dist; //Do the ride

   //update vehicle's position
   vh_pos[0] = rides[ride][2];
   vh_pos[1] = rides[ride][3];

   //If ride finished on time, add to the score
   if((*step) <= rides[ride][5]){
     int_score+=dist;
   }

   return int_score;
 }


static double getCarObjectiveValue(int car, int *n_rides, struct solution *s){
  struct problem *p = s->prob;
  int n = p->f + p->n;
  int i, step = 0, t = p->t;
  double score = .0, int_score;
  int vh_pos[2];

  vh_pos[0] = 0;
  vh_pos[1] = 0;
  i = s->cars[car]+1;

  for(; i<n && i<s->cars[car+1] && step < t && *n_rides > 0; i++){

    int_score = evaluateRide(&step, i, vh_pos, s);
    s->last_eval_ride = i;
    *n_rides--;
    score += int_score;
  }

  s->objv+=score;

  return score;
}

double *getObjectiveVector(double *objv, struct solution *s)
{
    /* solution is unfeasible, cannot evaluate it */
    double obj = 0.0;
    if (s->evalv) /* solution s is evaluated */
        *objv = s->objv;
    else { /* solution s is not evaluated */
      s->objv = 0;
      int n_rides = s->n_rides;
      // evaluate all solution cars
      int total_cars = s->prob->f;
      int i = 0;

      for(int vh = 0; vh<total_cars && n_rides>0; vh++){
        getCarObjectiveValue(vh, &n_rides, s);
      }

      *objv = -1*s->objv;
      s->evalv = 1;
    }

    return objv;
}

/*
 * Lower bound evaluation
 */

static double evaluateRideOptimistically(int i, struct solution *s){
 struct problem *p = s->prob;
 double int_score = 0.0;
 int **rides = p->rides;
 int dist;
 int ride = s->rides[i];

 //Give bonus if ride has started on time
 int_score+=p->b;

 // assume ride ends on time
 int_score += abs(rides[ride][0]-rides[ride][2]) + abs(rides[ride][1]-rides[ride][3]);

 return int_score;
}


double *getObjectiveLB(double *objLB, struct solution *s)
{
  struct problem *p = s->prob;
  int n = p->f + p->n;

  int vh_pos[2];
  vh_pos[0] = 0;
  vh_pos[1] = 0;

  double obj = 0.0;
  if (s->evalLB) /* solution s is evaluated */
      *objLB = -1*s->objLB;
  else { /* solution s is not evaluated */
      int n = s->prob->n + s->prob->f;
      getObjectiveVector(&obj, s);
      s->objLB = s->objv;

      int i = s->last_eval_ride+1;
      for(; i<n && s->rides[i] < s->prob->n; i++){
        s->objLB += evaluateRideOptimistically(i, s);
      }

      *objLB = -1*s->objLB;
      s->evalLB = 1;
  }
  return objLB;
}

/*
 * Modify a solution in place by applying a move to it
 */
struct solution *applyMove(struct solution *s, const struct move *v, const enum SubNeighbourhood nh)
{
    int i,j;
    int n = s->prob->n;
    int f = s->prob->f;
    switch (nh) {
    case ADD:
        // add to current car
        if(v->previousRide < n){
            // pos we will add to
            int pos = s->n_cars + s->n_rides;
            // find newRide's index
            for(j=pos; j<(n+f); j++){
                if(s->rides[j] == v->newRide)
                    break;
            }
            swap(s->rides, pos, j);
            s->n_rides++;
        } // add to new car
        else {
            // add car
            int pos = s->n_cars + s->n_rides;
            // find first car available to add's idnex
            for(j=pos; j<(n+f); j++){
                if(s->rides[j]>=n)
                    break;
            }
            swap(s->rides, pos, j);
            s->n_cars++;
            // add ride
            pos++;
            // find newRide's index
            for(j=pos; j<(n+f); j++){
                if(s->rides[j] == v->newRide)
                    break;
            }
            swap(s->rides, pos, j);
            s->n_rides++;
        }
        i = 0;
        break;
    case REMOVE:
        // remove from current car
        if(v->previousRide < n){
            s->n_rides--;
        } // remove car as well as ride
        // !! assuming car only had that one ride !!
        else {
            s->n_rides--;
            // find (from right to left) first position containing ride
            // so that end of s->rides contains only cars :)
            for(j=(n+f-1); j>=0; j--){
                if(s->rides[j]<n)
                    break;
            }
            // pos to add is n_cars+n_rides
            // we remove from the pos preceding it
            int pos = s->n_cars + s->n_rides - 1;
            swap(s->rides, pos, j);
            s->n_cars--;
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
    updateCarPositions(s);
    return s;
}

/*
 * Return true if a given solution is feasible or false if it is unfeasible
 */
int isFeasible(struct solution *s)
{
    return 1;
}

/*
 * Reset the enumeration of a given subneighbourhood of a solution
 */
struct solution *resetEnumMove(struct solution *s, const enum SubNeighbourhood nh)
{
    switch (nh) {
    case ADD:
        s->addMoveCounter = 0;
        break;
    case REMOVE:
        s->removeMoveCounter = 0;
        break;
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
    int a;
    switch (nh) {
    case ADD:
        a = (s->n_cars < s->prob->f && s->n_rides);
        return (a + 1)*(s->prob->n - s->n_rides);
        ->n_rides);
    case REMOVE:
        return 1;
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
    int aux;
    switch (st) {
    case PRESENT:
        if (!s->enumComponentState)
            return -1;
        aux = s->n_rides + s->n_cars - s->enumComponentState--;
        if (s->rides[aux] < s->prob->n) {
            return pairing(s->rides[aux - 1] < s->prob->n ? s->rides[aux - 1] : s->prob->n, s->rides[aux]);
        }
        else {
            s->enumComponentState--;
            return pairing(s->prob->n, s->rides[aux + 1]);
        }
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
        s->enumComponentState = s->n_rides + s->n_cars - 1;
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
    int n = s->prob->n;
    int f = s->prob->f;
    switch (nh) {
    case ADD:
        // all moves have been enumerated
        if(s->addMoveCounter >= getNeighbourhoodSize(s,nh))
            return NULL;
        // can't add anymore
        if(s->n_rides == n)
            return NULL;
        
        // first car must always have rides, so only N neighbourhood size
        if(s->n_rides==0){
            if(s->addMoveCounter >= n)
                return NULL;
            v->previousRide = -1; // previous ride was car
            v->newRide = s->rides[s->addMoveCounter+1];
        }// if all cars assigned, can only pool from remaining rides
        else if(s->n_cars==f){
            if(s->addMoveCounter >= n-s->n_rides)
                return NULL;
            // if s->ncars==f 
            // then we have already added the last car, and at least one ride (cf. applymove)
            // c r | r r r r 
            int pos = s->n_cars + s->n_rides;
            v->previousRide = s->rides[pos-1];
            v->newRide = s->rides[pos+s->addMoveCounter];
        }
        else {
            int pos = s->n_cars + s->n_rides;
            // add to current car
            if(s->addMoveCounter < n-s->n_rides){
                v->previousRide = s->rides[pos-1];
                v->newRide = s->rides[pos+s->addMoveCounter];
            } // add to new car
            else {
                v->previousRide = n;
                v->newRide = s->rides[pos+s->addMoveCounter];
            }
        }
        s->addMoveCounter++;
        break;
    case REMOVE:
        if(s->n_cars==1 && s->n_rides==0)
            return NULL;
        
        if(s->removeMoveCounter >= getNeighbourhoodSize(s,nh))
            return NULL;
        
        // if car has only one ride we remove the car as well
        // s->rides[pos-1] is car, s->rides[pos] is ride
        // s->rides[pos] can never be another car!
        // (since when we add car we also add ride, and first car always has ride)
        int pos = s->n_cars + s->n_rides - 1;
        if(s->rides[pos-1]>=n){ 
            // removal signaled with previousRide=n (not ride id)
            // if first car can't remove it! checked in applymove()
            v->previousRide = n;
            v->newRide = s->rides[pos];
        } // else remove car's last ride
        else {
            v->previousRide = s->rides[pos-1];
            v->newRide = s->rides[pos];
        }
        s->removeMoveCounter++;
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
    return v->componentId;
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
