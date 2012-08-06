#include "fixstring.h"

static const unsigned scales[6] = {
    1, 10, 100, 1000, 10000, 100000
};

void print_fix16_t(FILE *stream, fix16_t value, uint_fast8_t width, uint_fast8_t decimals)
{
    unsigned scale = scales[decimals];
    value = fix16_mul(value, scale);
    
    char fmt[5] = "%01d";
    fmt[2] = decimals + (value < 0 ? '2' : '1');
    char buf[20];
    int len = snprintf(buf, sizeof(buf), fmt, value);
    
    if (decimals) width--;
    width -= len;
    while (width-- > 0) fprintf(stream, " ");
    
    char tmp = buf[len - decimals];
    buf[len - decimals] = 0;
    
    fprintf(stream, "%s", buf);
    
    if (decimals != 0)
    {
        buf[len - decimals] = tmp;
        fprintf(stream, ".%s", &buf[len - decimals]);
    }
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
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ", ");
    print_fix16_t(stream, vector->x, 9, 4);
    fprintf(stream, ")");
}
