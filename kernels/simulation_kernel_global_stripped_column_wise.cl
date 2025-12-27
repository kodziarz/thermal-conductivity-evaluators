// #define DEBUG

// Cell types values (used in bufInBoard argument)
#define CELL_CONDUCTOR 0
#define CELL_GENERATOR 1
#define CELL_ADIABATIC 2
#define CELL_DRAIN 3

// kernel parameters set before compilation
// simulation parameters
#define STEPS_NUMBER 0
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

__kernel void simulate_heat(__global const int *bufInBoards,
                            __global float *maxTemperatures,
                            __global float *minTemperatures,
                            __global float *equilibriumMoments,
                            __global float *globalForegoingTemperatures,
                            __global float *globalNewTemperatures
#ifdef DEBUG
                            ,
                            __global float *debug
#endif
) {

  // #ifdef DEBUG
  //   debug[0] = -1.0; // started
  // #endif

  // calculate momory addresses
  int group_id = get_group_id(1);
  const int boardSize = WIDTH * HEIGHT;
  __global const int * currentBoard = bufInBoards + group_id * boardSize;
  __global float *foregoingTemperatures =
      globalForegoingTemperatures + group_id * boardSize;
  __global float *newTemperatures = globalNewTemperatures + group_id * boardSize;

  // calculate thread's coordinates
  const int stripsPerColumn = (HEIGHT - 2) / STRIP_LENGTH;
  const int col = get_local_id(1) + 1;
  const int stripIndex = get_local_id(0);
  const int global_id =
      (group_id * (WIDTH - 2) + col - 1) * stripsPerColumn + stripIndex;
  const int startRow = stripIndex * STRIP_LENGTH + 1;
  const int stripStartIndex = startRow * WIDTH + col;
  const int stripEndIndex = stripStartIndex + STRIP_LENGTH * WIDTH;

#ifdef DEBUG
  debug[global_id] = global_id; // ids calulated successfully
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
  debug[global_id] = -2.0; // alphas initiated
#endif

  // init temperatures tables

  for (int cellIndex = stripStartIndex; cellIndex < stripEndIndex;
       cellIndex += WIDTH) {

    int currentType = currentBoard[cellIndex];
    /*
    Branchless equivalent of:
    temperature = currentType == CELL_DRAIN
                                ? DRAIN_TEMPERATURE : START_TEMPERATURE;
    */
    {
      float temperature = (currentType == CELL_DRAIN) * DRAIN_TEMPERATURE +
                          (currentType != CELL_DRAIN) * START_TEMPERATURE;
      foregoingTemperatures[cellIndex] = temperature;
      newTemperatures[cellIndex] = temperature;
    }
    // other cell types are just border types, so they will not even use the
    // value
  }

  // initiate border - always with drain temperature (for adiabatic it does not
  // matter)
  if (stripIndex == 0) {
    foregoingTemperatures[stripStartIndex - WIDTH] = DRAIN_TEMPERATURE;
  } else if (stripIndex == stripsPerColumn - 1) {
    foregoingTemperatures[stripEndIndex + WIDTH] = DRAIN_TEMPERATURE;
  }

  if (col == 1) {
    for (int cellIndex = stripStartIndex - 1; cellIndex <= stripEndIndex - 1;
         cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = DRAIN_TEMPERATURE;
    }
  } else if (col == WIDTH - 2) {
    for (int cellIndex = stripStartIndex + 1; cellIndex <= stripEndIndex + 1;
         cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = DRAIN_TEMPERATURE;
    }
  }

  // corners are not initiated, since they are not used anyway

#ifdef DEBUG
  debug[global_id] = -10 - currentBoard[col + WIDTH]; // board copied
#endif

  barrier(CLK_GLOBAL_MEM_FENCE); // temperatures synchronization

  float maxT = foregoingTemperatures[stripStartIndex];
  float minT = foregoingTemperatures[stripStartIndex];
  int equilibriumMoment = 0;

  for (int step = 0; step < STEPS_NUMBER; step++) {

    // perform a step of a simulation

    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {

      const float foregoingT = foregoingTemperatures[cellIndex];
      float flow = 0;
      int currentColType = currentBoard[cellIndex];
      float beta = (currentColType == CELL_GENERATOR) * GENERATOR_BETA +
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

      float temperatureIncrease = DELTA_TIME * (flow + beta);

      newTemperatures[cellIndex] = foregoingT + temperatureIncrease;
    }

    // to that moment foregoingTemperatures have had proper foregoing values
    barrier(CLK_GLOBAL_MEM_FENCE);
    // from now on the foregoingTemperatures are actually undefined

    // copy data from newTemperatures to foregoingTemperatures
    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {
      float newTemperature = newTemperatures[cellIndex];
      foregoingTemperatures[cellIndex] = newTemperature;

      equilibriumMoment = (newTemperature > maxT) * step +
                          (newTemperature <= maxT) * equilibriumMoment;
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

  // find temperature minimum in equilibrium
  for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
       cellIndex += WIDTH) {
    float newTemperature = newTemperatures[cellIndex];
    /*
    Branchless equivalent of:
    if (newTemperature < minT)
      minT = newTemperature;
    */
    minT = (newTemperature < minT) * newTemperature +
           (newTemperature >= maxT) * minT;
  }

  // write to outputs of a strip
  maxTemperatures[global_id] = maxT;
  minTemperatures[global_id] = minT;
  equilibriumMoments[global_id] = equilibriumMoment;
}
