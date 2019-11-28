using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BloodFlow
{
    public struct Vector1
    {
        public Vector1(double _x)
        { x = _x;}
        public double x;

        public static Vector1 operator *(Vector1 v1, double d)
        {
            return new Vector1(v1.x * d);
        }

        public static Vector1 operator -(Vector1 v1, Vector1 v2)
        {
            return new Vector1(v1.x - v2.x);
        }

        public static Vector1 operator -(Vector1 v1)
        {
            return new Vector1(-v1.x);
        }

        public static bool operator <(Vector1 v1, Vector1 v2)
        {
            return v1.x < v2.x;
        }

        public static bool operator >(Vector1 v1, Vector1 v2)
        {
            return v1.x > v2.x;
        }

        public static Vector1 operator +(Vector1 v1, Vector1 v2)
        {
            return new Vector1(v1.x + v2.x);
        }

        static public double Distance(Vector1 v1, Vector1 v2)
        {
            return (double)Math.Abs(v1.x - v2.x);

        }

        static public double Dot(Vector1 v1, Vector1 v2)
        {
			return (v1.x * v2.x);
        }

        static public double Length(Vector1 v1)
        {
            return Math.Sqrt((v1.x * v1.x));
        }

        public double Normilize()
        {
            double l = (double)Math.Sqrt(x * x);
            if (l == 0)
            {
                x = 0;
                return 0;
            }
            x = x / l;
            return l;
        }

    };

    public delegate double SimpleFunction(double x);
    public delegate double GetBetaFunction(double R0, double elasticity);
    public delegate double MDFunction(double[] x);
}
