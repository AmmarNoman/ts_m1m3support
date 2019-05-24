
#define INTERVAL 6
#define VECTOR_SIZE 27
#define REPLICATE_ROW_SIZE 41254
#define REPLICATE_COL_SIZE 6
#define STRESS_SIZE REPLICATE_ROW_SIZE * VECTOR_SIZE * REPLICATE_COL_SIZE
#define ACTUATORS 156
#define BENDING_MODE_SIZE ACTUATORS * ACTUATORS

extern float coeff[VECTOR_SIZE];
extern float scale[VECTOR_SIZE];
extern unsigned int elemID[REPLICATE_ROW_SIZE];
extern unsigned int index[REPLICATE_ROW_SIZE];
extern float actuatorForces[ACTUATORS];
extern float bendingModesToActuatorForce[ACTUATORS*VECTOR_SIZE];
extern float actuatorForceToBendingModes[VECTOR_SIZE*ACTUATORS];
extern float stress[STRESS_SIZE];