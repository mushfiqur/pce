#ifndef _ENUMS_H
#define _ENUMS_H

// Shared enums used across classes

typedef enum {
	PCE,
	MONTE_CARLO
} SimType;

typedef enum{
	UNIFORM,
	GAUSSIAN
} RandVarDist;

typedef enum
{
	ADD,
	SUB,
	MULT,
	DIVIDE,
	INPUT_SIGNAL,
	INPUT_NOISE,
	DELAY,
	CONST,
	
	SINE_BLOCK,
	COSINE_BLOCK,
	FIR_BLOCK
} NodeType;

#endif