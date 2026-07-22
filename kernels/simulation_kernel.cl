enum Cell { CONDUCTOR = 0, GENERATOR = 1 };

enum BorderCell { ADIABATIC = 0, DRAIN = 1 };

__kernel void multiply_by_two(__global const uint8 *bufInBoards,
                              __global const uint8 *bufInBorders,
                              __global const double *tempreatureCalculations,
                              __global double *maxTemperatures,
                              const int height, const int width) {

  int global_id = get_global_id(0);
  int boardSize = height * width;
  int borderLength = 2 * (height + width) - 4;

  const uint8 *board = bufInBoards + global_id * boardSize;
  const uint8 *border = bufInBorders + global_id * borderLength;

  double result = all(board[0] == 1) * 100;

  maxTemperatures[global_id] = result;
  //   out[i] = in[i] * 2.0f;
}
