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

__kernel void simulate_heat(__global const int *bufInBoards, const int height,
                            const int width, __global float *maxTemperatures,
                            __local float *temperatures, __local int *board
#ifdef DEBUG
                            ,
                            __global float *debug
#endif
) {

#ifdef DEBUG
  debug[0] = -1.0; // ids calulated successfully
#endif

  // calculate addresses
  int group_id = get_group_id(1);
  const int boardSize = width * height;
  const int *currentBoard = bufInBoards + group_id * boardSize;

  int row = get_local_id(0);
  int col = get_local_id(1);
  int currentCellIndex = row * width + col;
  int global_id = group_id * boardSize + currentCellIndex;

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
  // initiate local memory
  board[currentCellIndex] = currentBoard[currentCellIndex];

  int currentCellType = board[currentCellIndex];
  double beta = GENERATOR_BETA;
  if (currentCellType == CELL_CONDUCTOR) {
    beta = CONDUCTOR_BETA;
  }
  // other cell types are just border types, so they will not even use the value

#ifdef DEBUG
  debug[global_id] = 10 + board[currentCellIndex]; // board copied
#endif

  temperatures[currentCellIndex] =
      currentCellType == CELL_DRAIN ? drainTemperature : startTemperature;

  work_group_barrier(CLK_LOCAL_MEM_FENCE);

  // if (col == 0 || col == width - 1 || row == 0 || row == height - 1) {
  //   // border point - nothing needs to be done

  //   for (int i = 0; i < SIMULATION_STEPS; i++) {
  //     work_group_barrier(CLK_LOCAL_MEM_FENCE); /*!< sychronization with end
  //     of
  //                        compute block (temperatures are not tampered with)*/
  //     // temperatures[currentCellIndex]
  //     work_group_barrier(CLK_LOCAL_MEM_FENCE); /*!<synchronization with end
  //     of
  //                                                 write-to-temperatures
  //                                                 block*/
  //   }
  //   return;
  // }

  float maxT = temperatures[currentCellIndex];
  bool isInterior =
      col != 0 && col != width - 1 && row != 0 && row != height - 1;
  for (int i = 0; i < numTimeSteps; i++) {
    float result = temperatures[currentCellIndex];
    if (isInterior) {
      // if not border cell
      const float foregoingT = temperatures[currentCellIndex];
      float flow = 0;

      int neighborIndex = currentCellIndex + width;
      flow += (temperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentCellType][board[neighborIndex]];

      neighborIndex = currentCellIndex - width;
      flow += (temperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentCellType][board[neighborIndex]];

      neighborIndex = currentCellIndex + 1;
      flow += (temperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentCellType][board[neighborIndex]];

      neighborIndex = currentCellIndex - 1;
      flow += (temperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentCellType][board[neighborIndex]];

      float temperatureIncrease = delta_time * (flow + beta);
      result = temperatures[currentCellIndex] + temperatureIncrease;
    }

    // to that moment temperatures proper foregoing values
    work_group_barrier(CLK_LOCAL_MEM_FENCE);
    // between this work_group_barriers temperatures have either foregoing or
    // next values (undefined)
    temperatures[currentCellIndex] = result;
    if (result > maxT)
      maxT = result;

    work_group_barrier(CLK_LOCAL_MEM_FENCE);
  }

  maxTemperatures[global_id] = maxT;
}
