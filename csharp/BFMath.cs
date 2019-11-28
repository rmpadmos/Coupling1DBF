using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BloodFlow
{
    public class myMath
    {
        public static double NewtonSolver(SimpleFunction func, double init_x, double eps, double dx)
        {
            double x = init_x;
			for (int i = 0; i < 10; i++)
            {
                double v0 = func(x);
                if (Math.Abs(v0) < eps)
                    return x;
                double dfdx = (func(x + dx) - v0) / dx;
				x = x - func(x) / dfdx;
            }
            return double.NaN;
        }

        public static double[] MDNewtonSolver(MDFunction[] funcs, bool[,] dep_matrix, double[] init_X, double eps, double[] dx)
        {
            int MAX_ITER_NUM = 10;
            int _info;
            alglib.densesolverreport _err_rep;
            
            int N = funcs.GetLength(0);
            double[,] m_jacobi = new double[N, N];
            double[] curr_x = (double[])init_X.Clone();
            double[] diff_x = (double[])init_X.Clone();
            double[] B = new double[N];

            int I = 0;
			for (I = 0; I < MAX_ITER_NUM; I++)
			{
				for (int i = 0; i < N; i++)
				{
					B[i] = funcs[i](curr_x);
					for (int j = 0; j < N; j++)
					{
						m_jacobi[i, j] = 0;
						if (dep_matrix[i, j])
						{
							curr_x[j] = curr_x[j] + dx[j];
                            m_jacobi[i, j] = funcs[i](curr_x);
                            m_jacobi[i, j] = (m_jacobi[i, j] - B[i]) / dx[j];
                            curr_x[j] = curr_x[j] - dx[j];
                        }
                    }
                }

                alglib.smp_rmatrixsolve(m_jacobi, N, B, out _info, out _err_rep, out diff_x);

                for (int i = 0; i < N; i++)
                {
                    curr_x[i] = curr_x[i] - diff_x[i];
                    diff_x[i] = curr_x[i];
                }

                double residual = 0;
                foreach (var f in funcs)
                    residual += Math.Abs(f(curr_x));



                if (residual < eps)
                    return curr_x;

                //      if (diff_x.Sum() == 0)
                //          for (int i = 0; i < diff_x.GetLength(0); i++)
                //              curr_x[i] = curr_x[i] * (1+1e-3);
            }

            return null;
        }
    }
}
