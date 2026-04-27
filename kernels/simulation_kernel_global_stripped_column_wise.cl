// #define DEBUG

// Cell types values (used in bufInBoard argument)
#define CELL_CONDUCTOR 0
#define CELL_GENERATOR 1
#define CELL_ADIABATIC 2
#define CELL_DRAIN 3

// kernel parameters set before compilation
// simulation parameters
#define STEPS_NUMBER 0
#define RESULTANT_POWER_TOL 0
#define GENERATOR_ALPHA 0
#define GENERATOR_BETA 0
#define CONDUCTOR_ALPHA 0
#define CONDUCTOR_BETA 0
#define DRAIN_ALPHA 0
#define DELTA_TIME 0
// board shape parameters
#define WIDTH 0
#define HEIGHT 0
#define THICKNESS 0
#define STRIP_LENGTH 1 // to avoid static code analysis errors due to division

typedef double simulation_value_t;
typedef int cell_type_t;
typedef int simulation_steps_index_t;

void initiate_temperatures(int depthStratumSize, int stripStartIndex,
                           int stripEndIndex, int col,
                           __global simulation_value_t *foregoingTemperatures,
                           __global simulation_value_t *newTemperatures,
                           __global const simulation_value_t *startTemperatures
#ifdef DEBUG
                           ,
                           int global_id, __global simulation_value_t *debug
#endif

) {
  // Border cells are not treated specially, because they need newTemperatures
  // initialization anyway - otherwise their final temperatures would be
  // misleading
  for (int depthStratumIndex = 0; depthStratumIndex < THICKNESS;
       depthStratumIndex++) {
    int depthStratumOffset = depthStratumIndex * depthStratumSize;
    // Take into account upper, lower, left-hand-side, right-hand-side, from
    // above and from below neighborhood rowCellIndex is the index of the cell
    // at the beginning of the row
    for (int rowCellIndex = depthStratumOffset + stripStartIndex - WIDTH - 1;
         rowCellIndex <= depthStratumOffset + stripEndIndex + WIDTH;
         rowCellIndex += WIDTH) {
      for (int cellIndex = rowCellIndex; cellIndex < rowCellIndex + 3;
           cellIndex++) {
        foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
        newTemperatures[cellIndex] = startTemperatures[cellIndex];
      }
    }
  }
}

void write_final_temperatures(int depthStratumSize, int stripStartIndex,
                              int stripEndIndex,
                              __global simulation_value_t *temperatures,
                              __global simulation_value_t *finalTemperatures

) {
  // Border cells are not treated specially, because they need newTemperatures
  // initialization anyway - otherwise their final temperatures would be
  // misleading
  for (int depthStratumIndex = 0; depthStratumIndex < THICKNESS;
       depthStratumIndex++) {
    int depthStratumOffset = depthStratumIndex * depthStratumSize;
    // Take into account upper, lower, left-hand-side, right-hand-side, from
    // above and from below neighborhood rowCellIndex is the index of the cell
    // at the beginning of the row
    for (int rowCellIndex = depthStratumOffset + stripStartIndex - WIDTH - 1;
         rowCellIndex <= depthStratumOffset + stripEndIndex + WIDTH;
         rowCellIndex += WIDTH) {
      for (int cellIndex = rowCellIndex; cellIndex < rowCellIndex + 3;
           cellIndex++) {
        finalTemperatures[cellIndex] = temperatures[cellIndex];
      }
    }
  }
}

__kernel void simulate_heat(
    // input data
    __global const cell_type_t *bufInBoards,
    __global const simulation_value_t *bufInStartTemperatures,
    // output data
    __global simulation_value_t *globalMaxTemperatures,          // per strip
    __global simulation_value_t *globalFinalTemperatures,        // per cell
    __global simulation_value_t *globalResultantPowers,          // per cell
    __global simulation_steps_index_t *globalEquilibriumMoments, // per strip
    // auxiliary tables for calculations
    __global simulation_value_t *globalForegoingTemperatures,
    __global simulation_value_t *globalNewTemperatures
#ifdef DEBUG
    ,
    __global simulation_value_t *debug
#endif
) {

  // #ifdef DEBUG
  //   // debug[0] = -1.0; // started
  //   printf("kernel run\n");
  // #endif

  // calculate memory addresses
  int group_id = get_group_id(1);
  const int depthStratumSize = WIDTH * HEIGHT;
  const int boardSize = depthStratumSize * THICKNESS;
  __global const cell_type_t *currentBoard = bufInBoards + group_id * boardSize;
  __global const simulation_value_t *startTemperatures =
      bufInStartTemperatures + group_id * boardSize;
  __global simulation_value_t *finalTemperatures =
      globalFinalTemperatures + group_id * boardSize;
  __global simulation_value_t *resultantPowers =
      globalResultantPowers + group_id * boardSize;
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
  printf("kernel %i run\n", global_id);
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
  printf("kernel %i initiated alphas\n", global_id);
#endif

  // init temperatures tables
  initiate_temperatures(depthStratumSize, stripStartIndex, stripEndIndex, col,
                        foregoingTemperatures, newTemperatures,
                        startTemperatures
#ifdef DEBUG
                        ,
                        global_id, debug
#endif
  );

#ifdef DEBUG
  debug[global_id] = -10 - currentBoard[col + WIDTH]; // board copied
  printf("kernel %i board copied\n", global_id);
#endif

  barrier(CLK_GLOBAL_MEM_FENCE); // temperatures synchronization

  simulation_value_t maxT = foregoingTemperatures[stripStartIndex];
  simulation_value_t minT = foregoingTemperatures[stripStartIndex];
  simulation_steps_index_t equilibriumMoment = 0;

  for (simulation_steps_index_t step = 0; step < STEPS_NUMBER; step++) {
    // perform a step of a simulation

#ifdef DEBUG
    debug[global_id] = -50 - global_id; // ids calulated successfully
    printf("kernel %i performing step %i\n", global_id, step);
#endif

    for (int depthStratumIndex = 1; depthStratumIndex < THICKNESS - 1;
         depthStratumIndex++) {

      int depthStratumOffset = depthStratumIndex * depthStratumSize;

      for (int cellIndex = depthStratumOffset + stripStartIndex;
           cellIndex <= depthStratumOffset + stripEndIndex;
           cellIndex += WIDTH) {

        const simulation_value_t foregoingT = foregoingTemperatures[cellIndex];
        simulation_value_t flow = 0;
        cell_type_t currentCellType = currentBoard[cellIndex];

        // power generated inside the cell
        /*
          Branchless equivalent of:
            simulation_value_t beta = currentCellType == CELL_GENERATOR ?
          GENERATOR_BETA : CONDUCTOR_BETA;
          */
        simulation_value_t beta =
            (currentCellType == CELL_GENERATOR) * GENERATOR_BETA +
            (currentCellType != CELL_GENERATOR) * CONDUCTOR_BETA;

        // order of neighbors: bottom, up, right, left, below, above

        int neighborIndex = cellIndex + WIDTH;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        neighborIndex = cellIndex - WIDTH;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        neighborIndex = cellIndex + 1;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        neighborIndex = cellIndex - 1;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        neighborIndex = cellIndex + depthStratumSize;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        neighborIndex = cellIndex - depthStratumSize;
        flow += (foregoingTemperatures[neighborIndex] - foregoingT) *
                ALPHAS[currentBoard[neighborIndex]][currentCellType];

        simulation_value_t resultantPower = flow + beta;
        resultantPowers[cellIndex] = resultantPower;

        simulation_value_t temperatureIncrease = DELTA_TIME * resultantPower;
        simulation_value_t newT = foregoingT + temperatureIncrease;
        newTemperatures[cellIndex] = foregoingT + temperatureIncrease;

        bool IsTolReached = resultantPower > RESULTANT_POWER_TOL;
        equilibriumMoment =
            IsTolReached * step + (!IsTolReached) * equilibriumMoment;
      }
    }

    // to that moment foregoingTemperatures have had proper foregoing values
    barrier(CLK_GLOBAL_MEM_FENCE);
    // from now on the foregoingTemperatures are actually undefined

    // copy data from newTemperatures to foregoingTemperatures
    for (int depthStratumIndex = 1; depthStratumIndex < THICKNESS - 1;
         depthStratumIndex++) {

      int depthStratumOffset = depthStratumIndex * depthStratumSize;

      for (int cellIndex = depthStratumOffset + stripStartIndex;
           cellIndex <= depthStratumOffset + stripEndIndex;
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
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
    // now foregoingTemperatures are well-defined again
  } // end of simulation

#ifdef DEBUG
  debug[global_id] = -2.0; // alphas initiated
  printf("kernel %i going to write final temperatures\n", global_id);
#endif

  // write to outputs of a strip
  // foregoingTemperatures are used on purpose for the edge-case of
  // STEPS_NUMBER == 0
  write_final_temperatures(depthStratumSize, stripStartIndex, stripEndIndex,
                           foregoingTemperatures, finalTemperatures);

  globalMaxTemperatures[global_id] = maxT;
  globalEquilibriumMoments[global_id] = equilibriumMoment;

#ifdef DEBUG
  debug[global_id] = -500; // ids calulated successfully
  printf("kernel %i finished successfully\n", global_id);
#endif
}
