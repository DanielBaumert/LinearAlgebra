# [Source](http://www.vilipetek.com/2014/12/13/polynomial-fitting-in-c/)

# Polynomial Fitting in C#
December 13, 2014 by vilipetek

![GitHub Logo](/img/PolyFitCSharpExample-300x200.png)

# Example

```csharp 
namespace DemoApplication
{
    class Program
    {
        const int DataSize = 1024;
        static void Main(string[] args)
        {
            
            double[] x = new double[DataSize];
            double[] y = new double[DataSize];
 
            // generate the data
            Random rand = new Random();
            for (int i = 0; i < DataSize; i++)
            {
                x[i] = i;
                y[i] = Math.Sin((double)i / 1024.0 * Math.PI * 2) + (rand.NextDouble() - 0.5);
            }
            
            // fit the data
            var polyfit = new PolyFit(x, y, 3);
            var fitted = polyfit.Fit(x);
        }
    }
}

```
Main class used for polynomial fitting is PolyFit. It relies on further internal classes but their listing is beyond the scope of this post.

```csharp
public class PolyFit
{
    /// <summary>
    /// Coefficients of a polynomial starting at the constant coefficient 
    /// and ending with the coefficient of power to the chosen order.
    /// </summary>
    public double[] Coeff { get; private set; }

    /// <summary>
    /// Finds the coefficients of a polynomial p(x) of degree n that fits the data, 
    /// p(x(i)) to y(i), in a least squares sense. The result p is a row vector of 
    /// length n+1 containing the polynomial coefficients in incremental powers.
    /// </summary>
    /// <param name="x">x axis values</param>
    /// <param name="y">y axis values</param>
    /// <param name="order">polynomial order including the constant</param>
    public PolyFit(double[] x, double[] y, int order)
    {
        // incrememnt the order to match matlab way
        double[,] matrixX = new double[x.Count(), ++order];
        double[,] matrixY = new double[x.Count(), 1];

        if (x.Length != y.Length)
        {
            throw new ArgumentException("x and y array lengths do not match!");
        }

        // copy y matrix
        for (int i = 0; i < y.Count(); i++)
        {
            matrixY[i, 0] = y[i];
        }

        // create the X matrix
        for (int row = 0; row < x.Count(); row++)
        {
            double nVal = 1.0f;
            for (int col = 0; col < order; col++)
            {
                matrixX[row, col] = nVal;
                nVal *= x[row];
            }
        }

        var matrixXt = matrixX.Transpose();
        var matrixXtX = matrixXt.Product(matrixX);
        var matrixXtY = matrixXt.Product(matrixY);

        var lu = new LUDecomposition(matrixXtX);
        Coeff = lu.Solve(matrixXtY).GetColumn(0).ToArray();
    }

    /// <summary>
    /// Calculates the value of a polynomial of degree n evaluated at x. The input argument 
    /// pCoeff is a vector of length n+1 whose elements are the coefficients in incremental 
    /// powers of the polynomial to be evaluated.
    /// </summary>
    /// <param name="x">Array of x values</param>
    /// <returns>Array of fitted y values</returns>
    public double[] Fit(double[] x)
    {
        double[] y = new double[x.Length];
        int pos = 0;

        foreach (double xval in x)
        {
            double xcoeff = 1;
            foreach (double coeffval in Coeff)
            {
                // multiply current x by a coefficient
                y[pos] += coeffval * xcoeff;
                // power up the X
                xcoeff *= xval;
            }
            pos++;
        }

        return y;
    }
}
```
