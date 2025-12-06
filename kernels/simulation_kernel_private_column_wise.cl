enum Cell { CONDUCTOR = 0, GENERATOR = 1 };

#define CELL_CONDUCTOR 0
#define CELL_GENERATOR 1
#define CELL_ADIABATIC 2
#define CELL_DRAIN 3
// OPTIMIZE global memory coalescence?

// #define DEBUG

const float deltaXY = 1.0;
const float deltaXY2 = deltaXY * deltaXY;
const int numTimeSteps = 10000; // Number of time steps to simulate

const float startTemperature = 280.0;
const float drainTemperature = 280.0;

// material params
// heat generator
const float g_m = 1000;
const float g_rho = 1000; // kg/m^3
const float g_k = 2000;
const float g_C_p = 50;  // J/kg*K
const float g_Q = 60000; // W

const float GENERATOR_ALPHA = g_k / (g_rho * g_C_p * deltaXY2);
const float GENERATOR_BETA = g_Q / (g_rho * g_C_p);

// heat conductor
const float c_m = 1000;
const float c_rho = 1000; // kg/m^3
const float c_k = 100000;
const float c_C_p = 5; // J/kg*K

const float CONDUCTOR_ALPHA = c_k / (c_rho * c_C_p * deltaXY2);
const float CONDUCTOR_BETA = 0;

const float d1 = g_m * g_C_p / (4 * g_k);
const float d2 = c_m * c_C_p / (4 * c_k);
const float delta_time = d1 < d2 ? d1 : d2;

const float DRAIN_ALPHA = 100;

#define STRIP_LENGTH 10
#define WIDTH 20
#define HEIGHT 40

// OPTIMIZE use struct as input
// typedef struct SimulationParams{

// }__attribute__((packed)) SimulationParams_t;

// void initialize_alphas(float ALPHAS[4][4], const float GENERATOR_ALPHA,
//                        const float CONDUCTOR_ALPHA, const float DRAIN_ALPHA)
//                        {
//   float c_g = (CONDUCTOR_ALPHA + GENERATOR_ALPHA) / 2;
//   float c_d = (CONDUCTOR_ALPHA + DRAIN_ALPHA) / 2;
//   float g_d = (GENERATOR_ALPHA + DRAIN_ALPHA) / 2;

//   for (int i = 0; i < 4; i++) {
//     ALPHAS[CELL_ADIABATIC][i] = 0;
//     ALPHAS[i][CELL_ADIABATIC] = 0;
//   }
//   ALPHAS[CELL_CONDUCTOR][CELL_GENERATOR] = c_g;
//   ALPHAS[CELL_GENERATOR][CELL_CONDUCTOR] = c_g;

//   ALPHAS[CELL_CONDUCTOR][CELL_DRAIN] = c_d;
//   ALPHAS[CELL_DRAIN][CELL_CONDUCTOR] = c_d;

//   ALPHAS[CELL_DRAIN][CELL_GENERATOR] = g_d;
//   ALPHAS[CELL_GENERATOR][CELL_DRAIN] = g_d;

//   ALPHAS[CELL_CONDUCTOR][CELL_CONDUCTOR] = CONDUCTOR_ALPHA;
//   ALPHAS[CELL_GENERATOR][CELL_GENERATOR] = GENERATOR_ALPHA;
//   ALPHAS[CELL_DRAIN][CELL_DRAIN] = 0;
// }

__kernel void simulate_heat(__global const int *bufInBoards,
                            __global float *maxTemperatures,
                            __local float *temperatures
#ifdef DEBUG
                            ,
                            __global float *debug
#endif
) {

#ifdef DEBUG
  debug[0] = -1.0; // ids calulated successfully
#endif

  // calculate addresses
  int group_id = get_group_id(0);
  const int boardSize = WIDTH * HEIGHT;
  const int *currentBoard = bufInBoards + group_id * boardSize;

  int col = get_local_id(0);
  int global_id = group_id * WIDTH + col;

#ifdef DEBUG
  debug[global_id] = 1.0; // ids calulated successfully
#endif

  // init ALPHAS
  float ALPHAS[4][4];
  // initialize_alphas(ALPHAS, GENERATOR_ALPHA, CONDUCTOR_ALPHA, DRAIN_ALPHA);
  float c_g = (CONDUCTOR_ALPHA + GENERATOR_ALPHA) / 2;
  float c_d = (CONDUCTOR_ALPHA + DRAIN_ALPHA) / 2;
  float g_d = (GENERATOR_ALPHA + DRAIN_ALPHA) / 2;

  for (int i = 0; i < 4; i++) {
    ALPHAS[CELL_ADIABATIC][i] = 0;
    ALPHAS[i][CELL_ADIABATIC] = 0;
  }
  ALPHAS[CELL_CONDUCTOR][CELL_GENERATOR] = c_g;
  ALPHAS[CELL_GENERATOR][CELL_CONDUCTOR] = c_g;

  ALPHAS[CELL_CONDUCTOR][CELL_DRAIN] = c_d;
  ALPHAS[CELL_DRAIN][CELL_CONDUCTOR] = c_d;

  ALPHAS[CELL_DRAIN][CELL_GENERATOR] = g_d;
  ALPHAS[CELL_GENERATOR][CELL_DRAIN] = g_d;

  ALPHAS[CELL_CONDUCTOR][CELL_CONDUCTOR] = CONDUCTOR_ALPHA;
  ALPHAS[CELL_GENERATOR][CELL_GENERATOR] = GENERATOR_ALPHA;
  ALPHAS[CELL_DRAIN][CELL_DRAIN] = 0;

#ifdef DEBUG
  debug[global_id] = 2.0; // alphas initiated
#endif

  bool isInterior = col != 0 && col != WIDTH - 1;

  // initiate local memory
  int leftColTypes[HEIGHT];
  int currentColTypes[HEIGHT];
  int rightColTypes[HEIGHT];

  float currentColBetas[HEIGHT];

  for (int cellIndex = col, row = 0; row < HEIGHT; cellIndex += WIDTH, row++) {

    int currentType = currentBoard[cellIndex];
    if (isInterior) {
      leftColTypes[row] = currentBoard[cellIndex - 1];
      currentColTypes[row] = currentType;
      rightColTypes[row] = currentBoard[cellIndex + 1];
      /*
      Branchless equivalent of:
      temperatures[cellIndex] = currentType == GENERATOR
                                 ? GENERATOR : CONDUCTOR_BETA;
      */
      currentColBetas[row] = (currentType == GENERATOR) * GENERATOR +
                             (currentType != GENERATOR) * CONDUCTOR_BETA;
    }
    /*
    Branchless equivalent of:
    temperatures[cellIndex] = currentType == CELL_DRAIN
                                ? drainTemperature : startTemperature;
    */
    temperatures[cellIndex] = (currentType == CELL_DRAIN) * drainTemperature +
                              (currentType != CELL_DRAIN) * startTemperature;
    // other cell types are just border types, so they will not even use the
    // value
  }

#ifdef DEBUG
  debug[global_id] = 10 + board[col]; // board copied
#endif

  work_group_barrier(CLK_LOCAL_MEM_FENCE); // temperatures synchronization

  float maxT = temperatures[col];

  for (int i = 0; i < numTimeSteps; i++) {

    // perform a step of a simulation
    float newColTemperatures[HEIGHT - 2];

    for (int cellIndex = WIDTH + col, row = 1; row < HEIGHT - 1;
         cellIndex += WIDTH, row++) {

      if (isInterior) {
        // if not border cell
        const float foregoingT = temperatures[cellIndex];
        float flow = 0;
        int currentColType = currentColTypes[row];

        int neighborIndex = cellIndex + WIDTH;
        flow += (temperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentColTypes[row + 1]][currentColType];

        neighborIndex = cellIndex - WIDTH;
        flow += (temperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentColTypes[row - 1]][currentColType];

        neighborIndex = cellIndex + 1;
        flow += (temperatures[neighborIndex] - foregoingT) *
                ALPHAS[rightColTypes[row]][currentColType];

        neighborIndex = cellIndex - 1;
        flow += (temperatures[neighborIndex] - foregoingT) *
                ALPHAS[leftColTypes[row - 1]][currentColType];

        float temperatureIncrease = delta_time * (flow + currentColBetas[row]);
        newColTemperatures[row - 1] = foregoingT + temperatureIncrease;
      }
    }

    // to that moment temperatures proper foregoing values
    work_group_barrier(CLK_LOCAL_MEM_FENCE);
    // temperatures values are undefined

    for (int cellIndex = col, row = 1; row < HEIGHT - 1;
         cellIndex += WIDTH, row++) {
      int newTemperature = newColTemperatures[row];
      temperatures[cellIndex] = newTemperature;
      /*
      Branchless equivalent of:
      if (newTemperature > maxT)
        maxT = newTemperature;
      */
      maxT = (newTemperature > maxT) * maxT +
             (newTemperature <= maxT) * newTemperature;
    }

    work_group_barrier(CLK_LOCAL_MEM_FENCE);
    // now temperatures are well-defined
  }

  maxTemperatures[global_id] = maxT;
}
