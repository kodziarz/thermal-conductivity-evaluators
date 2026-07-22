#pragma once

#include <vector>

#include "Types.h"
// enum Cell
// {
//     CONDUCTOR = 0,
//     GENERATOR = 1,
//     ADIABATIC = 2,
//     DRAIN = 3
// };

// using cell_type_t = int;

using SystemLayout_t = cell_type_t *;

namespace SystemLayout
{
    using namespace conductivity_evaluators;

    inline SystemLayout_t createGeneratorSystemLayout(int boardHeight, int boardWidth)
    {
        cell_type_t *result = new cell_type_t[boardHeight * boardWidth];

        // initialize with generator
        for (int i = 0; i < boardHeight * boardWidth; ++i)
            result[i] = Cell::GENERATOR;

        // top row and side columns -> adiabatic
        for (int row = 0; row < boardHeight - 1; row++)
        { // without last row
            for (int column = 0; column < boardWidth; column++)
            {
                if (row == 0 || column == 0 || column == boardWidth - 1)
                {
                    result[row * boardWidth + column] = Cell::ADIABATIC;
                }
            }
        }

        // bottom row: sides adiabatic, interior drains
        for (int column = 0; column < boardWidth; column++)
        {
            if (column == 0 || column == boardWidth - 1)
            {
                result[(boardHeight - 1) * boardWidth + column] = Cell::ADIABATIC;
            }
            else
            {
                result[(boardHeight - 1) * boardWidth + column] = Cell::DRAIN;
            }
        }

        return result;
    }

    inline SystemLayout_t createLeftConductorStripSystemLayout(int boardHeight, int boardWidth, int stripBredth = 1)
    {
        cell_type_t *result = new cell_type_t[boardHeight * boardWidth];

        // initialize with generator
        for (int i = 0; i < boardHeight * boardWidth; ++i)
            result[i] = Cell::GENERATOR;

        // top row and side columns -> adiabatic
        // left-most rows -> conductors
        for (int row = 0; row < boardHeight - 1; row++)
        { // without last row
            for (int column = 0; column < boardWidth; column++)
            {
                if (row == 0 || column == 0 || column == boardWidth - 1)
                {
                    result[row * boardWidth + column] = Cell::ADIABATIC;
                }
                else if (column <= stripBredth)
                {
                    result[row * boardWidth + column] = Cell::CONDUCTOR;
                }
            }
        }

        // bottom row: sides adiabatic, interior drains
        for (int column = 0; column < boardWidth; column++)
        {
            if (column == 0 || column == boardWidth - 1)
            {
                result[(boardHeight - 1) * boardWidth + column] = Cell::ADIABATIC;
            }
            else
            {
                result[(boardHeight - 1) * boardWidth + column] = Cell::DRAIN;
            }
        }

        return result;
    }

    template <typename Value_t>
    inline Value_t *rotateLeftSquareSystemLayoutBy90Deg(int boardLength, const Value_t *beingRotated)
    {
        Value_t *result = new Value_t[boardLength * boardLength];

        double translation = (boardLength - 1) / 2.0;

        for (int row = 0; row < boardLength; row++)
        {
            for (int column = 0; column < boardLength; column++)
            {
                double centeredRow = -(row - translation);
                double centeredColumn = column - translation;

                double centeredResultRow = centeredColumn;
                double centeredResultColumn = -centeredRow;

                int resultRow = -centeredResultRow + translation;
                int resultColumn = centeredResultColumn + translation;

                result[resultRow * boardLength + resultColumn] = beingRotated[row * boardLength + column];
            }
        }

        return result;
    }
}