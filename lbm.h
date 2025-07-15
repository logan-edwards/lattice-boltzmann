typedef struct vec2 {
	double x;
	double y;
} vec2;

typedef struct gridpoint {
	double weight[9];
	vec2 coordinate;
	double density;
	vec2 velocity;
	int boundary_condition;
} gridpoint;