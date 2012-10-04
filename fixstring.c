#include "fixstring.h"
#include <string.h>

void print_fix16_t(FILE *stream, fix16_t value, uint_fast8_t width, uint_fast8_t decimals)
{
    char buf[13];
    fix16_to_str(value, buf, decimals);
    
    uint_fast8_t len = strlen(buf);
    if (len < width)
    {
        width -= len;
        while (width-- > 0)
            fputc(' ', stream);
    }

    fputs(buf, stream);
}

void print_mf16(FILE *stream, const mf16 *matrix)
{
    if (matrix->errors)
    {
        fprintf(stream, "MATRIX ERRORS: %d\n", matrix->errors);
    }
    
    int row, column;
    for (row = 0; row < matrix->rows; row++)
    {
        for (column = 0; column < matrix->columns; column++)
        {
            fix16_t value = matrix->data[row][column];
            print_fix16_t(stream, value, 9, 4);
            fprintf(stream, " ");
        }
        fprintf(stream, "\n");
    }
}

void print_qf16(FILE *stream, const qf16 *quat)
{
    print_fix16_t(stream, quat->a, 9, 4);
    fprintf(stream, " ");
    print_fix16_t(stream, quat->b, 9, 4);
    fprintf(stream, "i ");
    print_fix16_t(stream, quat->c, 9, 4);
    fprintf(stream, "j ");
    print_fix16_t(stream, quat->d, 9, 4);
    fprintf(stream, "k");
}

void print_v3d(FILE *stream, const v3d *vector)
{
    fprintf(stream, "(");
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->y, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->z, 9, 4);
    fprintf(stream, ")");
}

void print_v2d(FILE *stream, const v2d *vector)
{
    fprintf(stream, "(");
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->y, 9, 4);
    fprintf(stream, ")");
}
