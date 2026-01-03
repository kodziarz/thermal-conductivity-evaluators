// #define DEBUG

// Cell types values (used in bufInBoard argument)
#define CELL_CONDUCTOR 0
#define CELL_GENERATOR 1
#define CELL_ADIABATIC 2
#define CELL_DRAIN 3

// kernel parameters set before compilation
// simulation parameters
#define STEPS_NUMBER 0
#define ETA 0
#define START_TEMPERATURE 0
#define DRAIN_TEMPERATURE 0
#define GENERATOR_ALPHA 0
#define GENERATOR_BETA 0
#define CONDUCTOR_ALPHA 0
#define CONDUCTOR_BETA 0
#define DRAIN_ALPHA 0
#define DELTA_TIME 0
// board shape parameters
#define WIDTH 0
#define HEIGHT 0
#define STRIP_LENGTH 1 // to avoid static code analysis errors due to division

typedef double simulation_value_t;
typedef int cell_type_t;
typedef int simulation_steps_index_t;

__kernel void
simulate_heat(__global const cell_type_t *bufInBoards,
              __global const simulation_value_t *bufInStartTemperatures,
              __global simulation_value_t *globalMaxTemperatures,
              __global simulation_value_t *globalFinalTemperatures,
              __global simulation_steps_index_t *globalEquilibriumMoments,
              __global simulation_value_t *globalForegoingTemperatures,
              __global simulation_value_t *globalNewTemperatures
#ifdef DEBUG
              ,
              __global simulation_value_t *debug
#endif
) {

  // #ifdef DEBUG
  //   debug[0] = -1.0; // started
  // #endif

  // calculate memory addresses
  int group_id = get_group_id(1);
  const int boardSize = WIDTH * HEIGHT;
  __global const cell_type_t *currentBoard = bufInBoards + group_id * boardSize;
  __global const simulation_value_t *startTemperatures =
      bufInStartTemperatures + group_id * boardSize;
  __global simulation_value_t *finalTemperatures =
      globalFinalTemperatures + group_id * boardSize;
  __global simulation_value_t *foregoingTemperatures =
      globalForegoingTemperatures + group_id * boardSize;
  __global simulation_value_t *newTemperatures =
      globalNewTemperatures + group_id * boardSize;

  // calculate thread's coordinates
  const int stripsPerColumn = (HEIGHT - 2) / STRIP_LENGTH;
  const int col = get_local_id(1) + 1;
  const int stripIndex = get_local_id(0);
  const int global_id = group_id * (WIDTH - 2) * stripsPerColumn +
                        stripIndex * (WIDTH - 2) + col - 1;
  const int startRow = stripIndex * STRIP_LENGTH + 1;
  const int stripStartIndex = startRow * WIDTH + col;
  const int stripEndIndex = stripStartIndex + (STRIP_LENGTH - 1) * WIDTH;

#ifdef DEBUG
  debug[global_id] = global_id; // ids calulated successfully
#endif

  // init ALPHAS
  simulation_value_t ALPHAS[4][4];
  // initialize_alphas(ALPHAS, GENERATOR_ALPHA, CONDUCTOR_ALPHA, DRAIN_ALPHA);
  simulation_value_t c_g = (CONDUCTOR_ALPHA + GENERATOR_ALPHA) / 2;
  simulation_value_t c_d = (CONDUCTOR_ALPHA + DRAIN_ALPHA) / 2;
  simulation_value_t g_d = (GENERATOR_ALPHA + DRAIN_ALPHA) / 2;

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
  debug[global_id] = -2.0; // alphas initiated
#endif

  // init temperatures tables

  for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
       cellIndex += WIDTH) {

    cell_type_t currentType = currentBoard[cellIndex];
    {
      foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
    }
  }

  // initiate border
  if (stripIndex == 0) {
    foregoingTemperatures[stripStartIndex - WIDTH] =
        startTemperatures[stripStartIndex - WIDTH];
    newTemperatures[stripStartIndex - WIDTH] =
        startTemperatures[stripStartIndex -
                          WIDTH]; // to read a meaningful border minimal
                                  // temperature later
  }
  if (stripIndex == stripsPerColumn - 1) {
    foregoingTemperatures[stripEndIndex + WIDTH] =
        startTemperatures[stripEndIndex + WIDTH];
    newTemperatures[stripStartIndex - WIDTH] =
        startTemperatures[stripStartIndex -
                          WIDTH]; // to read a meaningful border
                                  // minimal temperature later
  }

  if (col == 1) {
    for (int cellIndex = stripStartIndex - 1; cellIndex <= stripEndIndex - 1;
         cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
      newTemperatures[cellIndex] =
          startTemperatures[cellIndex]; // to read a meaningful border minimal
                                        // temperature later
    }
  }
  if (col == WIDTH - 2) {
    for (int cellIndex = stripStartIndex + 1; cellIndex <= stripEndIndex + 1;
         cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
      newTemperatures[cellIndex] =
          startTemperatures[cellIndex]; // to read a meaningful border minimal
                                        // temperature later
    }
  }

  // corners are not initiated, since they are not used anyway

#ifdef DEBUG
  debug[global_id] = -10 - currentBoard[col + WIDTH]; // board copied
#endif

  barrier(CLK_GLOBAL_MEM_FENCE); // temperatures synchronization

  simulation_value_t maxT = foregoingTemperatures[stripStartIndex];
  simulation_value_t minT = foregoingTemperatures[stripStartIndex];
  simulation_steps_index_t equilibriumMoment = 0;

  for (simulation_steps_index_t step = 0; step < STEPS_NUMBER; step++) {

    // perform a step of a simulation

    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {

      const simulation_value_t foregoingT = foregoingTemperatures[cellIndex];
      simulation_value_t flow = 0;
      cell_type_t currentColType = currentBoard[cellIndex];
      simulation_value_t beta =
          (currentColType == CELL_GENERATOR) * GENERATOR_BETA +
          (currentColType != CELL_GENERATOR) * CONDUCTOR_BETA;

      int neighborIndex = cellIndex + WIDTH;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentBoard[neighborIndex]][currentColType];

      neighborIndex = cellIndex - WIDTH;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentBoard[neighborIndex]][currentColType];

      neighborIndex = cellIndex + 1;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentBoard[neighborIndex]][currentColType];

      neighborIndex = cellIndex - 1;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
              ALPHAS[currentBoard[neighborIndex]][currentColType];

      simulation_value_t temperatureIncrease = DELTA_TIME * (flow + beta);
      simulation_value_t newT = foregoingT + temperatureIncrease;
      newTemperatures[cellIndex] = foregoingT + temperatureIncrease;

      simulation_value_t errorDenominator =
          (foregoingT > newT) * foregoingT + (foregoingT <= newT) * newT;
      bool didChange = errorDenominator != 0 &&
                       (temperatureIncrease / errorDenominator) >= ETA;
      equilibriumMoment = didChange * step + (!didChange) * equilibriumMoment;
    }

    // to that moment foregoingTemperatures have had proper foregoing values
    barrier(CLK_GLOBAL_MEM_FENCE);
    // from now on the foregoingTemperatures are actually undefined

    // copy data from newTemperatures to foregoingTemperatures
    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {

      simulation_value_t newTemperature = newTemperatures[cellIndex];
      foregoingTemperatures[cellIndex] = newTemperature;

      /*
            Branchless equivalent of:
            if (newTemperature > maxT)
              maxT = newTemperature;
            */
      maxT = (newTemperature > maxT) * newTemperature +
             (newTemperature <= maxT) * maxT;
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
    // now foregoingTemperatures are well-defined again
  } // end of simulation

  // write to outputs of a strip
  for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
       cellIndex += WIDTH) {
    finalTemperatures[cellIndex] = newTemperatures[cellIndex];
  }

  // write to outputs neigtbors if necessary
  if (stripIndex == 0) {
    finalTemperatures[stripStartIndex - WIDTH] =
        newTemperatures[stripStartIndex - WIDTH];
  }
  if (stripIndex == stripsPerColumn - 1) {
    finalTemperatures[stripEndIndex + WIDTH] =
        newTemperatures[stripEndIndex + WIDTH];
  }

  if (col == 1) {
    for (int cellIndex = stripStartIndex - 1; cellIndex <= stripEndIndex - 1;
         cellIndex += WIDTH) {
      finalTemperatures[cellIndex] = newTemperatures[cellIndex];
    }
  }
  if (col == WIDTH - 2) {
    for (int cellIndex = stripStartIndex + 1; cellIndex <= stripEndIndex + 1;
         cellIndex += WIDTH) {
      finalTemperatures[cellIndex] = newTemperatures[cellIndex];
    }
  }

  globalMaxTemperatures[global_id] = maxT;
  globalEquilibriumMoments[global_id] = equilibriumMoment;
}
