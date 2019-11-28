//#define FAST
#define SAFE

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

#if FAST
namespace BloodFlow
{
    unsafe public delegate double MDFunction_del(double* x);

    class delegateMDFunc : baseMDFunction
    {
        public delegateMDFunc(MDFunction_del _md_func_del)
        {
            md_func_del = _md_func_del;
        }

       unsafe public override double getVal(double* vec)
       {
           return md_func_del(vec);
       }
       private MDFunction_del md_func_del;
    }


    unsafe public delegate double _1DFunction_del(double x);
    class delegate1DFunc : base1DFunction
    {
        public delegate1DFunc(_1DFunction_del __1d_func_del)
        {
            _1d_func_del = __1d_func_del;
        }

        unsafe public override double getVal(double x)
        {
            return _1d_func_del(x);
        }
        private _1DFunction_del _1d_func_del;
    }
}
#endif 