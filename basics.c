#define PI 3.1415926535897932384626433832795
#define E 2.71828182845904523536028747135266245
#define Mu0 0.0000012566370614359172953850573533118

#define N 3e4    // 积分区间数
#define ESP 1e-6 // 第一类反常积分 IA与IB的被允许最大差值
#define CMAX 1e7 // 第一类反常积分 最大循环次数

long double integration(long double (*func)(long double), double a, double b, int n)
{
    // Simpson’s Three-Eighths Rule
    int count = 1;
    float f, fact2 = 0, fact1 = 0, integral, fact3 = 0;
    double h = (b - a) / n;
    double i = a + h;
    double f1 = func(a);
    double f2 = func(a + n * h);
    fact1 = f1 + f2;

    if ((n % 3) != 0)
    {
        // Invalid No. of Subintervals!
        return -1;
    }

    while (i <= b)
    {
        f = func(i);
        if (count % 3 == 0)
        {
            if (i <= n - 3)
                fact3 = fact3 + f;
        }
        else
        {
            if (i <= n - 1)
                fact2 = fact2 + f;
        }
        i = i + h;
        count++;
    }
    fact2 = fact2 * 3;
    fact3 = fact3 * 2;

    return (3 * h / 8) * (fact1 + fact2 + fact3);
}

long double integration_p2(long double (*func)(double, long double), double p, double a, double b, int n)
{
    // Simpson’s Three-Eighths Rule
    int count = 1;
    float f, fact2 = 0, fact1 = 0, integral, fact3 = 0;
    double h = (b - a) / n;
    double i = a + h;
    double f1 = func(p, a);
    double f2 = func(p, a + n * h);
    fact1 = f1 + f2;

    if ((n % 3) != 0)
    {
        // Invalid No. of Subintervals!
        return -1;
    }

    while (i <= b)
    {
        f = func(p, i);
        if (count % 3 == 0)
        {
            if (i <= n - 3)
                fact3 = fact3 + f;
        }
        else
        {
            if (i <= n - 1)
                fact2 = fact2 + f;
        }
        i = i + h;
        count++;
    }
    fact2 = fact2 * 3;
    fact3 = fact3 * 2;

    return (3 * h / 8) * (fact1 + fact2 + fact3);
}

long double integration_p3(long double (*func)(double, double, long double), double p1, double p2, double a, double b, int n)
{
    // Simpson’s Three-Eighths Rule
    int count = 1;
    float f, fact2 = 0, fact1 = 0, integral, fact3 = 0;
    double h = (b - a) / n;
    double i = a + h;
    double f1 = func(p1, p2, a);
    double f2 = func(p1, p2, a + n * h);
    fact1 = f1 + f2;

    if ((n % 3) != 0)
    {
        // Invalid No. of Subintervals!
        return -1;
    }

    while (i <= b)
    {
        f = func(p1, p2, i);
        if (count % 3 == 0)
        {
            if (i <= n - 3)
                fact3 = fact3 + f;
        }
        else
        {
            if (i <= n - 1)
                fact2 = fact2 + f;
        }
        i = i + h;
        count++;
    }
    fact2 = fact2 * 3;
    fact3 = fact3 * 2;

    return (3 * h / 8) * (fact1 + fact2 + fact3);
}

long double iEllipE(double k, long double x)
{
    return 1 / sqrt((1 - x * x) * (1 - k * k * x * x));
}
long double ellipE(long double k)
{
    return integration_p2(iEllipE, k, 0, 1, N);
}

long double iEllipK(double k, long double x)
{
    return sqrt((1 - k * k * x * x) / (1 - x * x));
}
long double ellipK(long double k)
{
    return integration_p2(iEllipK, k, 0, 1, N);
}

long double iJ0(double z, long double theta)
{
    return cos(z * cos(theta));
}
long double J0(long double z)
{
    return integration_p2(iJ0, z, 0, PI, N) / PI;
}

long double iJ1(double z, long double theta)
{
    return cos(theta - z * cos(theta));
}
long double J1(long double z)
{
    return integration_p2(iJ1, z, 0, PI, N) / PI;
}

long double iStruve0(double z, long double x)
{
    return sin(z * x) / sqrt(1 - x * x);
}
long double Struve0(long double z)
{
    return 2 * integration_p2(iStruve0, z, 0, 1, N) / PI;
}

long double iStruve1(double z, long double x)
{
    return sqrt(1 - x * x) * sin(z * x);
}
long double Struve1(long double z)
{
    return 2 * z * integration_p2(iStruve1, z, 0, 1, N) / PI;
}

inline long double U(double p, double x)
{
    return PI * (-J1(x) * Struve0(x) + p * J1(p * x) * Struve0(p * x) + J0(x) * Struve1(x) - p * J0(p * x) * Struve1(p * x)) / (2 * x * x);
}

long double iT(double p, double q, double x)
{
    return U(p, x) * U(p, x) * (q * x + pow(E, (-1 * q * x)) - 1);
}

long double calcL(double ri, double re, double l, double nc)
{
    unsigned long int upperLimit = 100;

    long double IA = 0;
    long double IB = 0;

    while (abs(IA - IB) <= ESP || upperLimit / 100 >= CMAX)
    {
        IA = IB;
        IB = 2 * PI * Mu0 * nc * nc * ri * ri * ri * ri * ri * integration_p3(iT, re / ri, l / ri, 0, upperLimit, N);

        upperLimit = upperLimit + 100;
    }

    return 0;
}

long double calcK(double Ra, double Rb, double d)
{
    return sqrt((4 * Ra * Rb) / ((Ra + Rb) * (Ra + Rb) + d * d));
}

long double calcM(double Ra, double Rb, double d)
{
    // TODO :@lru_cache()
    long double k = calcK(Ra, Rb, d);

    return Mu0 * sqrt(Ra * Rb) * ((2 / k - k) * ellipK(k) - (2 / k) * ellipE(k));
}

long double calcdM(double Ra, double Rb, double d)
{
    // TODO :@lru_cache()
    long double k = calcK(Ra, Rb, d);

    return (Mu0 * k * d * (2 * (1 -k * k) * ellipK(k) - (2 -k * k) * ellipE(k))) / (4 * (1 -k * k) * sqrt(Ra * Rb));
}

int main()
{
    return 0;
}