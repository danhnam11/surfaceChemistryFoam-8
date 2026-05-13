#include <stdio.h>

int readLastValue(const char* fileName, double* lastTime, double* lastValue)
{
    FILE* f = fopen(fileName, "r");
    char line[256];
    double t, v;
    int found = 0;

    if (!f)
    {
        printf("Cannot open %s\n", fileName);
        return 0;
    }

    while (fgets(line, sizeof(line), f))
    {
        if (line[0] == '#') continue;

        if (sscanf(line, "%lf %lf", &t, &v) == 2)
        {
            *lastTime = t;
            *lastValue = v;
            found = 1;
        }
    }

    fclose(f);
    return found;
}

int main()
{
    double tO2, tN2;
    double XO2out, XN2out;

    const double XO2in = 0.01;
    const double XN2in = 0.95;

    if (!readLastValue("postProcessing/outletXO2/0/surfaceFieldValue.dat", &tO2, &XO2out)) return 1;
    if (!readLastValue("postProcessing/outletXN2/0/surfaceFieldValue.dat", &tN2, &XN2out)) return 1;

    double conversion =
        1.0 - (XO2out/XN2out)/(XO2in/XN2in);

    printf("time_O2       = %.10g\n", tO2);
    printf("time_N2       = %.10g\n", tN2);
    printf("XO2out        = %.10e\n", XO2out);
    printf("XN2out        = %.10e\n", XN2out);
    printf("O2Conversion  = %.6f %%\n", conversion*100.0);

    return 0;
}
