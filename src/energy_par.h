/*
 *  energy_par.h
 *  The parameter value was obtained from http://www.cs.ubc.ca/labs/beta/Projects/RNA-Params 
 *  Reference: Mirela Andronescu et al. Computational approaches for RNA energy parameter estimation (RNA 2010)
 *  Created on: 2016/8/31
 *      Author: Tsukasa Fukunaga
 */

#ifndef ENERGY_PAR_H
#define ENERGY_PAR_H

#define GASCONST 1.98717  /* in [cal/K] */
#define K0  273.15
#define INF 1000000
#define TURN 3 
#define MAXLOOP 30

#include <string>

using namespace std;

static int temperature = 37;
static double kT = (temperature+K0)*GASCONST;
static double lxc37=107.856; /* parameter for logarithmic loop energy extrapolation*/

static int BP_pair[5][5]=
/* @  A  C  G  U*/
{{ 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5},
 { 0, 0, 0, 1, 0},
 { 0, 0, 2, 0, 3},
 { 0, 6, 0, 4, 0}};

/* rtype[pair[i][j]]:=pair[j][i] */
static int rtype[7] = {0, 2, 1, 4, 3, 6, 5}; 

static int hairpin37[31] = {
  INF, INF, INF, 364, 282, 297, 287, 259, 259, 272, 283,
  293, 303, 311, 319, 327, 334, 340, 346, 352, 358,
  363, 368, 373, 377, 382, 386, 390, 394, 398, 401};


static int mismatchH37[7][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -90,  -35,  -13,    3,   13}, /* A@  AA  AC  AG  AU */
   { -90,  -12,  -17,  -98,  -22}, /* C@  CA  CC  CG  CU */
   { -90,  -51,  -24,   -1,    3}, /* G@  GA  GC  GG  GU */
   { -90,  -16,  -57,  -97,  -38}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { -70,    1,  -22,  -65,  110}, /* A@  AA  AC  AG  AU */
   { -70,   11,  -27,  -17,    0}, /* C@  CA  CC  CG  CU */
   { -70,  -67,    1,  -45,   36}, /* G@  GA  GC  GG  GU */
   { -70,    8,  -39,  -55,  -92}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   60,   75,   43,  180}, /* A@  AA  AC  AG  AU */
   {   0,   50,   44,   70,   -1}, /* C@  CA  CC  CG  CU */
   {   0,  -92,   15,   39,  114}, /* G@  GA  GC  GG  GU */
   {   0,   25,    7,   67,   36}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   43,   52,   87,  109}, /* A@  AA  AC  AG  AU */
   {   0,   52,   35,  -10,   40}, /* C@  CA  CC  CG  CU */
   {   0,  -26,  -28,    5,   92}, /* G@  GA  GC  GG  GU */
   {   0,   42,  -38,  -14,  -30}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   41,   77,   65,  117}, /* A@  AA  AC  AG  AU */
   {   0,   -2,    8,   44,   43}, /* C@  CA  CC  CG  CU */
   {   0,  -31,   23,   15,   60}, /* G@  GA  GC  GG  GU */
   {   0,   29,   -6,   32,  -10}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   17,   52,   76,   68}, /* A@  AA  AC  AG  AU */
   {   0,   53,   36,   37,   11}, /* C@  CA  CC  CG  CU */
   {   0,  -39,   46,   16,   68}, /* G@  GA  GC  GG  GU */
   {   0,   50,  -16,   46,   29}}/* U@  UA  UC  UG  UU */
};

static int stack37[7][7] =
/*          CG     GC     GU     UG     AU     UA  */
{ {  INF,   INF,   INF,   INF,   INF,   INF,   INF},
  {  INF,  -133,  -207,  -146,   -37,  -139,  -132},
  {  INF,  -207,  -205,  -150,   -91,  -129,  -123},
  {  INF,  -146,  -150,   -22,   -67,   -81,   -58},
  {  INF,   -37,   -91,   -67,   -38,   -14,    -2},
  {  INF,  -139,  -129,   -81,   -14,   -84,   -69},
  {  INF,  -132,  -123,   -58,    -2,   -69,   -68}};

static int bulge37[31] = {
  INF, 281, 157, 201, 288, 298, 272, 288, 303, 315,
  327, 337, 346, 355, 363, 370, 377, 384, 390, 396,
  401, 407, 412, 416, 421, 425, 430, 434, 438, 442,
  445};

static int TerminalAU = 56;

static int internal_loop37[31] = {
  INF, INF, INF, INF,  83, 118,  90, 106, 121, 133,
  145, 155, 164, 173, 181, 188, 195, 202, 208, 214,
  219, 225, 230, 234, 239, 243, 248, 252, 256, 260,
  263};

static int mismatchI37[7][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0,  -50,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0,  -50,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -45}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,    0,    0,  -50,    0}, /* A@  AA  AC  AG  AU */
   {   0,    0,    0,    0,    0}, /* C@  CA  CC  CG  CU */
   {   0,  -50,    0,    0,    0}, /* G@  GA  GC  GG  GU */
   {   0,    0,    0,    0,  -45}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   63,   63,   13,   63}, /* A@  AA  AC  AG  AU */
   {   0,   63,   63,   63,   63}, /* C@  CA  CC  CG  CU */
   {   0,   13,   63,   63,   63}, /* G@  GA  GC  GG  GU */
   {   0,   63,   63,   63,   18}},/* U@  UA  UC  UG  UU */
  { /* UG */
  {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   63,   63,   13,   63}, /* A@  AA  AC  AG  AU */
   {   0,   63,   63,   63,   63}, /* C@  CA  CC  CG  CU */
   {   0,   13,   63,   63,   63}, /* G@  GA  GC  GG  GU */
   {   0,   63,   63,   63,   18}},/* U@  UA  UC  UG  UU */
  { /* AU */
  {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   63,   63,   13,   63}, /* A@  AA  AC  AG  AU */
   {   0,   63,   63,   63,   63}, /* C@  CA  CC  CG  CU */
   {   0,   13,   63,   63,   63}, /* G@  GA  GC  GG  GU */
   {   0,   63,   63,   63,   18}},/* U@  UA  UC  UG  UU */
  { /* UA */
  {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   {   0,   63,   63,   13,   63}, /* A@  AA  AC  AG  AU */
   {   0,   63,   63,   63,   63}, /* C@  CA  CC  CG  CU */
   {   0,   13,   63,   63,   63}, /* G@  GA  GC  GG  GU */
   {   0,   63,   63,   63,   18}},/* U@  UA  UC  UG  UU */
};

static int ML_closing37 = 315;
static int ML_intern37 = 15;
static int ML_BASE37 = -2;

static int dangle5_37[8][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF}, /* no pair */
   { INF,   -8,    0,    0,    0}, /* CG  (stacks on C) */
   { INF,  -10,    0,    0,    0}, /* GC  (stacks on G) */
   { INF,    0,    0,    0,    0}, /* GU */
   { INF,    0,    0,    0,    0}, /* UG */
   { INF,  -10,    0,    0,    0}, /* AU */
   { INF,   -3,    0,    0,    0}, /* UA */
   {   0,    0,     0,    0,   0}  /*  @ */
};

/* 3' dangling ends (unpaired base stacks on second paired base */
static int dangle3_37[8][5]=
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   { INF,  -10,   -9,  -51,  -14},  /* CG  (stacks on G) */
   { INF,  -41,    0,  -46,  -35},  /* GC */
   { INF,  -10,    0,  -50,   -2},  /* GU */
   { INF,  -10,  -40,  -60,    0},  /* UG */
   { INF,  -10,  -13,  -25,   -7},  /* AU */
   { INF,  -10,  -30,  -44,   -7},  /* UA */
   {   0,    0,     0,    0,   0}   /*  @ */
};

static int MAX_NINIO = 300;
static int F_ninio37 = 50;

#endif
