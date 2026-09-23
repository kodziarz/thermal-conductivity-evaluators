#pragma once

#include <vector>

#include "Types.h"

using SystemLayout_t = int *;

namespace SystemLayout
{
    using namespace conductivity_evaluators;

    enum CellType
    {
        CONDUCTOR = 0,
        GENERATOR = 1,
        ADIABATIC = 2,
        DRAIN = 3
    };

    struct MaterialProperties
    {
        simulation_value_t k;
        simulation_value_t invC;
        simulation_value_t qGen;
    };

    inline void cellTypeLayoutToProperties(
        const int *layout, int size,
        const MaterialProperties &conductorProps,
        const MaterialProperties &generatorProps,
        const MaterialProperties &adiabaticProps,
        const MaterialProperties &drainProps,
        simulation_value_t *k_out,
        simulation_value_t *invC_out,
        simulation_value_t *qGen_out)
    {
        for (int i = 0; i < size; i++)
        {
            MaterialProperties props;
            switch (layout[i])
            {
            case CONDUCTOR:
                props = conductorProps;
                break;
            case GENERATOR:
                props = generatorProps;
                break;
            case ADIABATIC:
                props = adiabaticProps;
                break;
            case DRAIN:
                props = drainProps;
                break;
            default:
                props = {0, 0, 0};
                break;
            }
            k_out[i] = props.k;
            invC_out[i] = props.invC;
            qGen_out[i] = props.qGen;
        }
    }

    inline SystemLayout_t createGeneratorSystemLayout(int boardHeight, int boardWidth)
    {
        int *result = new int[boardHeight * boardWidth];

        // initialize with generator
        for (int i = 0; i < boardHeight * boardWidth; ++i)
            result[i] = GENERATOR;

        // top row and side columns -> adiabatic
        for (int row = 0; row < boardHeight - 1; row++)
        { // without last row
            for (int column = 0; column < boardWidth; column++)
            {
                if (row == 0 || column == 0 || column == boardWidth - 1)
                {
                    result[row * boardWidth + column] = ADIABATIC;
                }
            }
        }

        // bottom row: sides adiabatic, interior drains
        for (int column = 0; column < boardWidth; column++)
        {
            if (column == 0 || column == boardWidth - 1)
            {
                result[(boardHeight - 1) * boardWidth + column] = ADIABATIC;
            }
            else
            {
                result[(boardHeight - 1) * boardWidth + column] = DRAIN;
            }
        }

        return result;
    }

    inline SystemLayout_t createLeftConductorStripSystemLayout(int boardHeight, int boardWidth, int stripBredth = 1)
    {
        int *result = new int[boardHeight * boardWidth];

        // initialize with generator
        for (int i = 0; i < boardHeight * boardWidth; ++i)
            result[i] = GENERATOR;

        // top row and side columns -> adiabatic
        // left-most rows -> conductors
        for (int row = 0; row < boardHeight - 1; row++)
        {
            for (int column = 0; column < boardWidth; column++)
            {
                if (row == 0 || column == 0 || column == boardWidth - 1)
                {
                    result[row * boardWidth + column] = ADIABATIC;
                }
                else if (column <= stripBredth)
                {
                    result[row * boardWidth + column] = CONDUCTOR;
                }
            }
        }

        // bottom row: sides adiabatic, interior drains
        for (int column = 0; column < boardWidth; column++)
        {
            if (column == 0 || column == boardWidth - 1)
            {
                result[(boardHeight - 1) * boardWidth + column] = ADIABATIC;
            }
            else
            {
                result[(boardHeight - 1) * boardWidth + column] = DRAIN;
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
