typedef struct vec2 {
	double x;
	double y;
} vec2;

typedef struct gridpoint {
	vec2 coordinate;
	vec2 velocity;
	double density;
	double f[9];
	int has_boundary_condition;	// might be unnecessary, handling can occur in func loop
} gridpoint;

typedef struct color {
	unsigned int r;
	unsigned int g;
	unsigned int b;
	unsigned int alpha;
} color;