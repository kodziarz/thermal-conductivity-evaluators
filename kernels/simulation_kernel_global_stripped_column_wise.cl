// #define DEBUG

// kernel parameters set before compilation
// simulation parameters
#define STEPS_NUMBER 0
#define ETA 0
#define DELTA_TIME 0
#define DRAIN_TEMPERATURE 0
// board shape parameters
#define WIDTH 0
#define HEIGHT 0
#define STRIP_LENGTH 1 // to avoid static code analysis errors due to division

typedef double simulation_value_t;
typedef int simulation_steps_index_t;

__kernel void
simulate_heat(__global const simulation_value_t *bufInStartTemperatures,
              __global simulation_value_t *globalMaxTemperatures,
              __global simulation_value_t *globalFinalTemperatures,
              __global simulation_steps_index_t *globalEquilibriumMoments,
              __global simulation_value_t *globalForegoingTemperatures,
              __global simulation_value_t *globalNewTemperatures,
              __global const simulation_value_t *bufK,
              __global const simulation_value_t *bufInvC,
              __global const simulation_value_t *bufQGen
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
  __global const simulation_value_t *startTemperatures =
      bufInStartTemperatures + group_id * boardSize;
  __global simulation_value_t *finalTemperatures =
      globalFinalTemperatures + group_id * boardSize;
  __global simulation_value_t *foregoingTemperatures =
      globalForegoingTemperatures + group_id * boardSize;
  __global simulation_value_t *newTemperatures =
      globalNewTemperatures + group_id * boardSize;
  __global const simulation_value_t *cellK =
      bufK + group_id * boardSize;
  __global const simulation_value_t *cellInvC =
      bufInvC + group_id * boardSize;
  __global const simulation_value_t *cellQGen =
      bufQGen + group_id * boardSize;

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

  // init temperatures tables

  for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
       cellIndex += WIDTH) {
    foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
  }

  // initiate border
  if (stripIndex == 0) {
    foregoingTemperatures[stripStartIndex - WIDTH] =
        startTemperatures[stripStartIndex - WIDTH];
    newTemperatures[stripStartIndex - WIDTH] =
        startTemperatures[stripStartIndex -
                          WIDTH]; // for correctness of returned border T
  }
  if (stripIndex == stripsPerColumn - 1) {
    foregoingTemperatures[stripEndIndex + WIDTH] =
        startTemperatures[stripEndIndex + WIDTH];
    newTemperatures[stripEndIndex + WIDTH] =
        startTemperatures[stripEndIndex +
                          WIDTH]; // for correctness of returned border T
  }

  if (col == 1) {
    // starting from row above and ending a row below, to cover corners
    // (theoretically not necessary, but done for the sake of correctness of
    // returned border minimal temperature)
    for (int cellIndex = stripStartIndex - 1 - WIDTH;
         cellIndex <= stripEndIndex + WIDTH - 1; cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
      newTemperatures[cellIndex] =
          startTemperatures[cellIndex]; // for correctness of returned border T
    }
  }
  if (col == WIDTH - 2) {
    // starting from row above and ending a row below, to cover corners
    // (theoretically not necessary, but done for the sake of correctness of
    // returned border minimal temperature)
    for (int cellIndex = stripStartIndex - WIDTH + 1;
         cellIndex <= stripEndIndex + WIDTH + 1; cellIndex += WIDTH) {
      foregoingTemperatures[cellIndex] = startTemperatures[cellIndex];
      newTemperatures[cellIndex] =
          startTemperatures[cellIndex]; // for correctness of returned border T
    }
  }

  // corners are not initiated, since they are not used anyway

#ifdef DEBUG
  debug[global_id] = -10.0; // board copied
#endif

  barrier(CLK_GLOBAL_MEM_FENCE); // temperatures synchronization

  simulation_value_t maxT = foregoingTemperatures[stripStartIndex];
  simulation_steps_index_t equilibriumMoment = 0;

  for (simulation_steps_index_t step = 0; step < STEPS_NUMBER; step++) {
    // perform a step of a simulation

    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {

      const simulation_value_t foregoingT = foregoingTemperatures[cellIndex];
      simulation_value_t flow = 0;

      simulation_value_t alpha_current = cellK[cellIndex] * cellInvC[cellIndex];
      simulation_value_t beta = cellQGen[cellIndex] * cellInvC[cellIndex];

      simulation_value_t alpha_neighbor;
      simulation_value_t mutual_alpha;

      int neighborIndex = cellIndex + WIDTH;
      alpha_neighbor = cellK[neighborIndex] * cellInvC[neighborIndex];
      mutual_alpha = (cellK[cellIndex] == 0 || cellK[neighborIndex] == 0)
                         ? 0
                         : (alpha_current + alpha_neighbor) / 2;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) * mutual_alpha;

      neighborIndex = cellIndex - WIDTH;
      alpha_neighbor = cellK[neighborIndex] * cellInvC[neighborIndex];
      mutual_alpha = (cellK[cellIndex] == 0 || cellK[neighborIndex] == 0)
                         ? 0
                         : (alpha_current + alpha_neighbor) / 2;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) * mutual_alpha;

      neighborIndex = cellIndex + 1;
      alpha_neighbor = cellK[neighborIndex] * cellInvC[neighborIndex];
      mutual_alpha = (cellK[cellIndex] == 0 || cellK[neighborIndex] == 0)
                         ? 0
                         : (alpha_current + alpha_neighbor) / 2;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) * mutual_alpha;

      neighborIndex = cellIndex - 1;
      alpha_neighbor = cellK[neighborIndex] * cellInvC[neighborIndex];
      mutual_alpha = (cellK[cellIndex] == 0 || cellK[neighborIndex] == 0)
                         ? 0
                         : (alpha_current + alpha_neighbor) / 2;
      flow += (foregoingTemperatures[neighborIndex] - foregoingT) * mutual_alpha;

      simulation_value_t temperatureIncrease = DELTA_TIME * (flow + beta);
      simulation_value_t newT = foregoingT + temperatureIncrease;
      newTemperatures[cellIndex] = foregoingT + temperatureIncrease;

      simulation_value_t absForegoingT = fabs(foregoingT);
      simulation_value_t absNewT = fabs(newT);
      simulation_value_t errorDenominator =
          (absForegoingT > absNewT) * absForegoingT +
          (absForegoingT <= absNewT) * absNewT;
      bool didChange = errorDenominator != 0 &&
                       fabs(temperatureIncrease / errorDenominator) >= ETA;
      equilibriumMoment = didChange * step + (!didChange) * equilibriumMoment;
#ifdef DEBUG
      if (step == 49000) {
        debug[global_id] = ETA; // board copied
      }
#endif
    }

    // to that moment foregoingTemperatures have had proper foregoing values
    barrier(CLK_GLOBAL_MEM_FENCE);
    // from now on the foregoingTemperatures are actually undefined

    // copy data from newTemperatures to foregoingTemperatures
    // and reset DRAIN cells (invC == 0) to drainTemperature
    for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
         cellIndex += WIDTH) {

      simulation_value_t newTemperature = newTemperatures[cellIndex];

      // DRAIN cells: invC == 0, force to drainTemperature
      simulation_value_t isDrain = (cellInvC[cellIndex] == 0);
      newTemperature = isDrain * DRAIN_TEMPERATURE + (!isDrain) * newTemperature;

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
  // foregoingTemperatures are used on purpose for the edge-case of
  // STEPS_NUMBER == 0
  for (int cellIndex = stripStartIndex; cellIndex <= stripEndIndex;
       cellIndex += WIDTH) {
    finalTemperatures[cellIndex] = foregoingTemperatures[cellIndex];
  }

  // write to outputs neighbors if necessary
  if (stripIndex == 0) {
    finalTemperatures[stripStartIndex - WIDTH] =
        foregoingTemperatures[stripStartIndex - WIDTH];
  }
  if (stripIndex == stripsPerColumn - 1) {
    finalTemperatures[stripEndIndex + WIDTH] =
        foregoingTemperatures[stripEndIndex + WIDTH];
  }

  if (col == 1) {
    for (int cellIndex = stripStartIndex - WIDTH - 1;
         cellIndex <= stripEndIndex + WIDTH - 1; cellIndex += WIDTH) {
      finalTemperatures[cellIndex] = foregoingTemperatures[cellIndex];
    }
  }
  if (col == WIDTH - 2) {
    for (int cellIndex = stripStartIndex - WIDTH + 1;
         cellIndex <= stripEndIndex + WIDTH + 1; cellIndex += WIDTH) {
      finalTemperatures[cellIndex] = foregoingTemperatures[cellIndex];
    }
  }

  globalMaxTemperatures[global_id] = maxT;
  globalEquilibriumMoments[global_id] = equilibriumMoment;
}
