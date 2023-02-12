/**  
 * A recreation of the `java.lang.Math` class.
 * This does not use any inbuilt or external methods from other classes. Use according to your wishes.
 * @author SrcyDev
 * @version 1.0.1
 * 
 * Know Issues: pow(),exp(),tanh(),cosh(),sinh() and some others  i.e. hyperbolic,logarithmic and exponential functions do not work properly.
 *
 * Note - This class's method may not be as accurate than the inbuilt `java.lang.Math' class. Also, it may provide incorrect results.
 */

 public class MathFunctions {
    public static final double PI = 3.141592653589793;
    public static final double E = 2.718281828459045;
    public static final double SQRT2 = 1.414213562373095;
    public static final double SQRT3 = 1.732050807568877;
    public static final double LOG2E = 1.4426950408889634;
    public static final double LOG10E = 0.4342944819032518;
    public static final double LN2 = 0.6931471805599453;
    public static final double LN10 = 2.302585092994046;
    public static final int MAX_BYTE = 127;
    public static final int MIN_BYTE = -128;
    public static final int MAX_SHORT = 32767;
    public static final int MIN_SHORT = -32768;
    public static final int MAX_INT = 2147483647;
    public static final int MIN_INT = -2147483648;
    public static final long MAX_LONG = 9223372036854775807L;
    public static final long MIN_LONG = -9223372036854775808L;
    public static final float MAX_FLOAT = 3.4028235E38f;
    public static final float MIN_FLOAT = 1.4E-45f;
    public static final double MAX_DOUBLE = 1.7976931348623157E308d;
    public static final double MIN_DOUBLE = 4.9E-324d;
    public static final char MAX_CHAR = '\uFFFF';
    public static final char MIN_CHAR = '\u0000';
    public static final boolean TRUE = true;
    public static final boolean FALSE = false;

    public static int abs(int x) {
        return (x < 0) ? -x : x;
    }

    public static long abs(long x) {
        return (x < 0) ? -x : x;
    }

    public static float abs(float x) {
        return (x < 0) ? -x : x;
    }

    public static double abs(double x) {
        return (x < 0) ? -x : x;
    }

    public static double max(double x, double y) {
        return x > y ? x : y;
    }

    public static double min(double x, double y) {
        return x < y ? x : y;
    }

    public static double round(double x) {
        double r = x % 1;

        if (r >= 0.5) {
            x += (1 - r);
        } else {
            x -= r;
        }

        return x;
    }

    public static double random(long seed) {
        seed = (seed * 25214903917L + 11) & ((1L << 48) - 1);
        return (seed >> 12) / (double)(1L << 48);
    }

    public static double sqrt(double x) {
        if (x < 0) {
            throw new IllegalArgumentException("Cannot compute square root of negative number");
        }

        double result = x;
        double epsilon = 1e-15;
        while (abs(result - x / result) > epsilon * result) {
            result = (x / result + result) / 2.0;
        }
        return result;
    }

    public static double cbrt(double x) {
        if (x == 0) {
            return 0;
        }

        double guess = x / 3;
        double lastGuess;
        do {
            lastGuess = guess;
            guess = (2 * lastGuess + x / (lastGuess * lastGuess)) / 3;
        } while (abs(guess - lastGuess) > 1e-15);

        return guess;
    }

    public static double fractionalRoot(double x, double y) {
        // Use the Newton-Raphson method to iteratively improve the estimate of the root
        double guess = y;
        double lastGuess = 0;
        double tolerance = 1e-10; // Set the tolerance for the error between successive guesses
        while (abs(guess - lastGuess) > tolerance) {
            lastGuess = guess;
            guess = ((y - 1) * guess + x / pow(guess, y - 1)) / y;
        }
        return guess;
    }

    public static double pow(double base, double exponent) {
        return exp(exponent * log10(base));
    }

    public static long factorial(int n) {
        long result = 1;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }

    public static double signum(double x) {
        if (Double.isNaN(x)) {
          return Double.NaN;
        } else if (x == Double.POSITIVE_INFINITY) {
          return 1;
        } else if (x == Double.NEGATIVE_INFINITY) {
          return -1;
        } else if (x < 0) {
          return -1;
        } else if (x > 0) {
          return 1;
        } else {
          return 0;
        }
    }

    public static double ceil(double x) {
        if (x > (int) x) {
          return (int) x + 1;
        }
        return (int) x;
    }

    public static double copySign(double x, double y) {
        long xBits = Double.doubleToLongBits(x);
        long yBits = Double.doubleToLongBits(y);
        yBits = (yBits & 0x7FFFFFFFFFFFFFFFL) | (xBits & 0x8000000000000000L);
        return Double.longBitsToDouble(yBits);
    }

    public static double nextAfter(double start, double direction) {
        long bits = Double.doubleToLongBits(start);
        long sign = bits & (1L << 63);
        long exponent = bits & ((1L << 52) - 1);
        long mantissa = bits & ((1L << 52) - 1);
    
        if (exponent == (1L << 52) - 1) {
            // Handle infinite and NaN values
            if (mantissa == 0) {
                // Return infinity with the same sign as start
                return Double.longBitsToDouble(sign | (1L << 52) - 1);
            } else {
                // Return NaN
                return Double.longBitsToDouble((1L << 52) - 1);
            }
        }
    
        if (direction > start) {
            // Return the next larger value
            if (mantissa == (1L << 52) - 1) {
                // The mantissa is already at the maximum value, so we need to carry the one to the exponent
                exponent += 1;
                mantissa = 0;
            } else {
                // Increment the mantissa
                mantissa += 1;
            }
        } else {
            // Return the next smaller value
            if (mantissa == 0) {
                // The mantissa is already at the minimum value, so we need to borrow the one from the exponent
                exponent -= 1;
                mantissa = (1L << 52) - 1;
            } else {
                // Decrement the mantissa
                mantissa -= 1;
            }
        }
    
        // Reassemble the bits and return the result
        return Double.longBitsToDouble(sign | exponent | mantissa);
    }

    public static double nextUp(double x) {
        if (x == Double.POSITIVE_INFINITY) {
            return x;
        }
        if (x == Double.NEGATIVE_INFINITY) {
            return Double.longBitsToDouble(1L << 63);
        }
        if (x == 0.0) {
            return MIN_DOUBLE;
        }
        if (x == -0.0) {
            return -MIN_DOUBLE;
        }
        if (x < 0) {
            return Double.longBitsToDouble(Double.doubleToRawLongBits(x) - 1);
        }
        return Double.longBitsToDouble(Double.doubleToRawLongBits(x) + 1);
    }

    public static double nextDown(double x) {
        if (x == Double.NEGATIVE_INFINITY) {
            return x;
        }
        if (x == Double.POSITIVE_INFINITY) {
            return Double.longBitsToDouble((1L << 63) - 1);
        }
        if (x == 0.0) {
            return -Double.longBitsToDouble((1L << 63) - 1);
        }
        if (x == -0.0) {
            return Double.longBitsToDouble((1L << 63) - 1);
        }
        if (x > 0) {
            return Double.longBitsToDouble(Double.doubleToRawLongBits(x) - 1);
        }
        return Double.longBitsToDouble(Double.doubleToRawLongBits(x) + 1);
    }

    public static int floor(double value) {
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            throw new IllegalArgumentException("Cannot floor NaN or infinite values");
        }
        if (value < Integer.MIN_VALUE) {
            return Integer.MIN_VALUE;
        }
        if (value > Integer.MAX_VALUE) {
            return Integer.MAX_VALUE;
        }
        return (int) value;
    }

    public static int floorDiv(int x, int y) {
        int result = x / y;
        if (x % y < 0) {
            result--;
        }
        return result;
    }

    public static int rint(double value) {
        return (int) round(value);
    }

    public static double hypot(double x, double y) {
        double xAbs = abs(x);
        double yAbs = abs(y);
        if (xAbs > yAbs) {
          double temp = yAbs / xAbs;
          return xAbs * sqrt(1 + temp * temp);
        } else if (yAbs > 0) {
          double temp = xAbs / yAbs;
          return yAbs * sqrt(1 + temp * temp);
        } else {
          return 0;
        }
    }

    public static double ulp(double d) {
        if (d == 0.0) {
            return MIN_DOUBLE;
        }
        return abs(d - nextUp(d));
    }

    public static int getExponent(double num) {
        int exponent = 0;
        while (num >= 10) {
            num /= 10;
            exponent++;
        }
        return exponent;
    }

    public static double IEEEremainder(double dividend, double divisor) {
        double remainder = dividend - (floor(dividend / divisor) * divisor);
        return remainder;
    }

    public static int addExact(int x, int y) {
        int sum = x + y;
        if (x > 0 && y > 0 && sum < 0) {
            throw new ArithmeticException("Integer overflow");
        }
        if (x < 0 && y < 0 && sum > 0) {
            throw new ArithmeticException("Integer underflow");
        }
        return sum;
    }

    public static int subtractExact(int x, int y) {
        int result = x - y;
        if (x < 0 && y > 0 && result > x) {
            throw new ArithmeticException("Integer overflow");
        }
        if (x > 0 && y < 0 && result < x) {
            throw new ArithmeticException("Integer overflow");
        }
        return result;
    }

    public static int multiplyExact(int x, int y) {
        int result = 0;
        while (y != 0) {
            if (y % 2 != 0) {
                result += x;
            }
            x *= 2;
            y /= 2;
        }
        return result;
    }

    public static int incrementExact(int number) {
        int result = number + 1;
        if (result < number) {
            throw new ArithmeticException("Integer overflow");
        }
        return result;
    }

    public static int decrementExact(int number) {
        return number - 1;
    }

    public static int negateExact(int x) {
        int neg = 0;
        int d = (x > 0) ? -1 : 1;
        while (x != 0) {
            neg += d;
            x += d;
        }
        return neg;
    }

    public static int toIntExact(long value) {
        if ((value & 0xFFFFFFFF00000000L) == 0) {
            return (int) value;
        }
        if (value > 0) {
            return Integer.MAX_VALUE;
        } else {
            return Integer.MIN_VALUE;
        }
    }

    public static double log(double base, double x) {
        // Initialize the sum to 0 and the power to x - 1
        double sum = 0;
        double power = x - 1;
        // Set the tolerance for the error between successive terms
        double tolerance = 1e-10;
        // Initialize the counter to 1
        int n = 1;
        // Calculate the sum of the series until the error between successive terms is less than the tolerance
        while (abs(power) > tolerance) {
            sum += power / n;
            power *= x - 1;
            n++;
        }
        // Divide the sum by the logarithm of the base to get the logarithm of x to the base base
        return sum / log(base, base);
    }

    public static double log10(double x) {
        if (Double.isNaN(x) || Double.isInfinite(x)) {
          throw new IllegalArgumentException("Number must be real and finite");
        }
        
        if (x <= 0) {
            throw new IllegalArgumentException("Number must be positive");
        }
        // Initialize the sum to 0 and the power to x - 1
        double sum = 0;
        double power = x - 1;
        // Set the tolerance for the error between successive terms
        double tolerance = 1e-15;
        // Initialize the counter to 1
        int n = 1;
        // Calculate the sum of the series until the error between successive terms is less than the tolerance
        while (abs(power) > tolerance) {
            sum += power / n;
            power *= x - 1;
            n++;
        }
        // Divide the sum by the natural logarithm of 10 to get the logarithm of x to base 10
        return sum / 2.303;
    }

    public static double exp(double x) {
        double result = 1.0;
        for (int i = 1; i < 20; i++) {
            result += pow(x, i) / factorial(i);
        }
        return result;
    }

    public static double expm1(double x) {
        double result = 0;
        double term = 1;
        int i = 1;
        while (term > 0.00001 || term < -0.00001) {
            result += term;
            term *= x / i;
            i++;
        }
        return result;
    }

    public static double log1p(double x) {
        double result = 0;
        double term = 1;
        int i = 1;
        while (term > 1e-15) {
            term *= (x / i);
            result += term;
            i++;
        }
        return result;
    }

    public static double toDegrees(double radians) {
        return radians * 180 / PI;
    }
  
      public static double toRadians(double degrees) {
        return degrees * PI / 180;
    }

    public static double sin(double x) {
        double sin = 0;
        if (x < 0) {
            x = -x;
            sin = -sin(x);
        } else if (x > Math.PI / 2) {
            x = Math.PI - x;
            sin = sin(x);
        } else {
            double x2 = x * x;
            sin = x * (1 + x2 * (-1.0 / 6 + x2 * (1.0 / 120 - x2 / 5040)));
        }
        return sin;
    }

  
    public static double cos(double x) {
        double g = sin(x);
        return sqrt(1 - (g * g)); // Pythagorean Identity : sin ^ 2(x) + cos ^ 2 (x) = 1
    }
  
    public static double tan(double x)
    {
        double g = sin(x);
        return g / sqrt(1- (g * g)); // As tan(x) = sin(x) / cos(x)
    }
  
      public static double asin(double x) {
        if (x < -1 || x > 1) {
            throw new IllegalArgumentException("Invalid input for asin function: " + x);
        }
          
        return 1/ sin(x);
    }
  
    public static double acos(double x)
    {
        if (x < -1 || x > 1) {
          throw new IllegalArgumentException("Invalid input for acos function: " + x);
        }
  
        return 1/cos(x);
  
    }
  
    public static double atan(double x)
    {
        if (x < -1 || x > 1) {
          throw new IllegalArgumentException("Invalid input for acos function: " + x);
        }
  
        return 1/tan(x);
    }
  
      public static double atan2(double y, double x) {
        // Handle special cases
        if (x == 0) {
            if (y > 0) return MathFunctions.PI / 2;
            if (y == 0) return 0;
            return -MathFunctions.PI / 2;
        }
  
        // Calculate the value of atan(y/x) using a Taylor series expansion
        //double absX = abs(x);
        //double absY = abs(y);
        double z = y / x;
        double z2 = z * z;
        double result = z;
        for (int i = 1; i < 10; i++) {
            result += ((i % 2 == 1) ? -1 : 1) * z2 / (2 * i + 1);
            z2 *= z;
        }
  
        // Adjust the result based on the signs of x and y
        if (x < 0) {
            if (y < 0) return result - MathFunctions.PI;
            return result + MathFunctions.PI;
        } else {
            return result;
        }
    }

    public static double sinh(double x) {
        return (exp(x) - exp(-x)) / 2;
    }
  
    public static double cosh(double x) {
        return (exp(x) + exp(-x)) / 2;
    }
  
    public static double tanh(double x) {
        return sinh(x) / cosh(x);
    }
}
