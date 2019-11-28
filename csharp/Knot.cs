#define SAFE
//#define FAST

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;



namespace BloodFlow
{
    public class Knot
    {
        public Knot(VascularNode _core_node, double start_time)
        {
            core_node = _core_node;
            nodes = (VascularNode[])core_node.neighbours.ToArray().Clone();

            int L = nodes.GetLength(0);
            
            velocity = new double[L];
            pressure = new double[L];
            lumen_sq = new double[L];

            for (int i = 0; i < core_node.neighbours.Count; i++)
            {               
                lumen_sq[i] = core_node.neighbours[i].lumen_sq_0;
            }

            DefineSigns();

            current_time = start_time;
            previous_time = current_time - 1e-3;
        }

        public void DefineSigns()
        {
            int L = nodes.GetLength(0);
            Vector1[] dir_vector1 = new Vector1[L];
            Vector1[] dir_vector2 = new Vector1[L];
            Vector1[] dir_vector3 = new Vector1[L];
            Vector1[] dir_vector4 = new Vector1[L];
            v_sign = new int[L];
            v_sign_1 = new int[L];

            for (int i = 0; i < L; i++)
            {
                //dir_vector1[i] = core_node.position - nodes[i].position;// change this to constant? remember sign
                if (Math.Abs(nodes[i].position.x) < 1e-10) //redefined direction
                    dir_vector1[i] = new Vector1(-1.0);
                else
                    dir_vector1[i] = new Vector1(1.0);
                dir_vector2[i] = nodes[i].dir_vector;
                dir_vector3[i] = nodes[i].position - nodes[i].neighbours.Last().position;
                dir_vector4[i] = nodes[i].neighbours.Last().dir_vector;
            }

            for (int i = 0; i < L; i++)
            {
                v_sign[i] = Math.Sign(Vector1.Dot(dir_vector1[i], dir_vector2[i]));
                v_sign_1[i] = Math.Sign(Vector1.Dot(dir_vector3[i], dir_vector4[i]))*v_sign[i];
            }
        }

        public int NaNControl()
        {
            foreach (var n in nodes)
                if (double.IsNaN(n.pressure) || double.IsNaN(n.velocity))
                    return core_node.id;
            return -1;
        }

        public virtual void reset()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                velocity[i] = 0;
                pressure[i] = 0;
                lumen_sq[i] = nodes[i].lumen_sq_0;
                nodes[i].velocity = velocity[i];
                nodes[i].pressure = pressure[i];
                nodes[i].lumen_sq = nodes[i].lumen_sq_0;
            } 
        }

        public virtual void doCoupling(double dt) { }

        public virtual void holdChtr(double dt) { }

        public VascularNode[] nodes;
        public VascularNode core_node;


        public int[] v_sign;
        public int[] v_sign_1;
        public double[] velocity;
        public double[] pressure;
        public double[] lumen_sq;
        public double[] lumen_sq_0;

        public double current_time;
        public double previous_time;
    }

    public class
           StandartKnot : Knot
    {
        public StandartKnot(Knot _knot, GetBetaFunction getElasticBeta)
            : base(_knot.core_node, _knot.current_time)
        {
            int L = nodes.GetLength(0);
            chrt_func = new MDFunction[L];
            energy_conservation_func = new MDFunction[L - 1];
            funcs = new MDFunction[2 * L];

            wall_thickhess = new double[L];
            lumen_sq_0 = new double[L];

            beta_1 = new double[L];
            chrt_b = new double[L];
            chrt_f = new double[L];
            c_dst  = new double[L];
            dep_matrix = new bool[2 * L, 2 * L];
            prev_velocity = new double[L];
            g_energy = new double[L];

#if FAST
            nl_system = new NewtonSolver(2 * L);
#endif


            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    dep_matrix[i, j] = false;

            for (int i = 0; i < L; i++)
            {
                double R0 = Math.Sqrt(nodes[i].lumen_sq_0 / Math.PI);
                beta_1[i] = getElasticBeta(R0, nodes[i].elasticity) / nodes[i].lumen_sq_0;
                wall_thickhess[i] = GlobalDefs.getBoileauWallThickness(R0, nodes[i].elasticity);
                lumen_sq_0[i] = nodes[i].lumen_sq_0;
                prev_velocity[i] = nodes[i].velocity;
                chrt_b[i] = 0;//-(4 * Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f /  GlobalDefs.BLOOD_DENSITY));
                chrt_f[i] = 0;// (4 * Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f /  GlobalDefs.BLOOD_DENSITY));               
                c_dst[i] = Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f / GlobalDefs.BLOOD_DENSITY);
                g_energy[i] = 0;//GlobalDefs.BLOOD_DENSITY * Vector3.Dot(GlobalDefs.DOWN, nodes[i].position - GlobalDefs.ZERO_POINT) * GlobalDefs.GRAVITY;
            }

            for (int i = 0; i < L; i++)
                pressure[i] = GlobalDefs.DIASTOLIC_PRESSURE;



            int count = 0;
            unsafe
            {
                for (int i = 0; i < L; i++)
                {
                    int I = i;
#if FAST
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v > 0)
                            return Math.Abs(v) + 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_b[I];
                    };

                    baseMDFunction f1 = new delegateMDFunc(f1_del);
                    nl_system.addFunc(f1);

                    nl_system.setDetMatrixEl(count, 2 * I, true);
                    nl_system.setDetMatrixEl(count, 2 * I + 1, true);
#endif

                    chrt_func[i] = delegate(double[] args) //v1,l1; v2,l2 ...
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v > 0)
                            return Math.Abs(v) + 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_b[I];
                    };

                    
                    funcs[count] = chrt_func[i];

                    dep_matrix[count, 2 * I] = true;
                    dep_matrix[count, 2 * I + 1] = true;


                    count++;
                }
            }

            unsafe
            {
#if FAST
                MDFunction_del f1_del = delegate(double* args)
                {
                    double summ_flux = 0;
                    for (int i = 0; i < L; i++)
                        summ_flux += args[0 + i * 2] * args[1 + i * 2];

                    return summ_flux;
                };
                baseMDFunction f1 = new delegateMDFunc(f1_del);
                nl_system.addFunc(f1);

                for (int i = 0; i < 2 * L; i++)              
                    nl_system.setDetMatrixEl(count, i, true);
                
#endif

                mass_conservation_func = delegate(double[] args)
                {
                    double summ_flux = 0;
                    for (int i = 0; i < L; i++)
                    {
                        double v = args[0 + i * 2];
                        double l = args[1 + i * 2];

                        summ_flux += v * l;
                    }
                    return summ_flux;
                };


                funcs[count] = mass_conservation_func;
                for (int i = 0; i < 2 * L; i++)               
                    dep_matrix[count, i] = true;
               
               
            };

            count++;

            unsafe
            {
                for (int i = 1; i < L; i++)
                {
                    int I = i;
#if FAST
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v0 = args[0];
                        double p0 = beta_1[0] * (Math.Sqrt(args[1]) - Math.Sqrt(nodes[0].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;

                        double v = args[0 + I * 2];
                        double p = beta_1[I] * (Math.Sqrt(args[1 + 2 * I]) - Math.Sqrt(nodes[I].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;
                        return GlobalDefs.BLOOD_DENSITY * (v0 * v0 - v * v) / 2 + p0 - p + g_energy[0] - g_energy[I];
                    };
                    baseMDFunction f1 = new delegateMDFunc(f1_del);
                    nl_system.addFunc(f1);

                    nl_system.setDetMatrixEl(count, 0, true);
                    nl_system.setDetMatrixEl(count, 1, true);
                    nl_system.setDetMatrixEl(count, 2 * I, true);
                    nl_system.setDetMatrixEl(count, 2 * I + 1, true);
#endif

                    energy_conservation_func[i - 1] = delegate(double[] args)
                    {
                        double v0 = args[0];
                        double p0 = beta_1[0] * (Math.Sqrt(args[1]) - Math.Sqrt(nodes[0].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;

                        double v = args[0 + I * 2];
                        double p = beta_1[I] * (Math.Sqrt(args[1 + 2 * I]) - Math.Sqrt(nodes[I].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;
						//return 0;
						return GlobalDefs.BLOOD_DENSITY * (v0 * v0 - v * v) / 2 + p0 - p + g_energy[I];

                    };
                    
                    funcs[count] = energy_conservation_func[I - 1];

                    dep_matrix[count, 0] = true;
                    dep_matrix[count, 1] = true;
                    dep_matrix[count, 2 * I] = true;
                    dep_matrix[count, 2 * I + 1] = true;

                    count++;
                }

#if FAST
                us_init_X = (double*)Marshal.AllocHGlobal(2 * L * sizeof(double));
                us_solution = (double*)Marshal.AllocHGlobal(2 * L * sizeof(double));

                for (int i = 0; i < 2 * L; i += 2)
                {
                    us_init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                    us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                    nl_system.setDxVectorEl(i, 1e-12f);
                    nl_system.setDxVectorEl(i + 1, 1e-12f);
                }
#endif
            }         

            dX = new double[2 * L];
            for (int i = 0; i < 2 * L; i += 2)
            {
                dX[i] = 1e-12f;
                dX[i + 1] = 1e-12f;
            }

        }

        public override void reset()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                velocity[i] = 0;
                pressure[i] = GlobalDefs.DIASTOLIC_PRESSURE;
                lumen_sq[i] = nodes[i].lumen_sq_0;
                nodes[i].velocity = velocity[i];
                nodes[i].pressure = pressure[i];
                nodes[i].lumen_sq = nodes[i].lumen_sq_0;
            }

            unsafe
            {
                for (int i = 0; i < 2 * L; i += 2)
                {
                    us_init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                    us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                }
            }
        }


        unsafe public override void doCoupling(double dt)
        {
            previous_time = current_time;
            current_time = current_time + dt;

            int L = nodes.GetLength(0);

            //if(core_node.id==359)
                //L = nodes.GetLength(0);


            for (int i = 0; i < L; i++)
                nodes[i].lumen_sq = this.calcLumen_sq(i, nodes[i].pressure);

#if SAFE            
            double[] solution = new double[2 * L];
            double[] init_X   = new double[2 * L];     
            for(int i=0; i<2*L; i+=2)
            {
                init_X[i  ] = nodes[i/2].velocity*v_sign[i/2];
                init_X[i+1] = nodes[i/2].lumen_sq;
                double wave_speed = 4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[i/2]);
                // Comment of nex two rows gives new algo for knot couplong
                chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + wave_speed;//4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[i/2]);
                chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - wave_speed;//4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[i/2]); 
            }

            solution = myMath.MDNewtonSolver(funcs, dep_matrix, init_X, 1e-6, dX);

            double av_pressure = 0;
            double av_flux_in = 0;
            double av_flux_out = 0;
            double av_lumen_in = 0;
            double av_lumen_out = 0;

            for (int i = 0; i < 2 * L; i += 2)
            {
                nodes[i / 2].velocity = solution[i] * v_sign[i / 2];
                nodes[i / 2].lumen_sq = solution[i + 1];
                nodes[i / 2].pressure = beta_1[i / 2] * (Math.Sqrt(nodes[i / 2].lumen_sq) - Math.Sqrt(nodes[i / 2].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;

                av_pressure += nodes[i / 2].pressure;
                    if (solution[i] >= 0)
                    {
                        av_flux_in += Math.Abs(solution[i + 1] * solution[i]);
                        av_lumen_in += solution[i + 1];
                    }
                    else
                    {
                        av_flux_out += Math.Abs(solution[i + 1] * solution[i]);
                        av_lumen_out += solution[i + 1];
                    }
                core_node.pressure = av_pressure / L;
                core_node.lumen_sq = av_lumen_in;
                core_node.velocity = av_flux_in / av_lumen_in;
            }

            core_node.velocity = core_node.neighbours.Last().velocity;
            core_node.lumen_sq = core_node.neighbours.Last().lumen_sq;
#endif


#if FAST
            unsafe
            {
                //double* us_init_X   = (double*)Marshal.AllocHGlobal(2* L * sizeof(double));
                //double* us_solution = (double*)Marshal.AllocHGlobal(2* L * sizeof(double));  


                for (int i = 0; i < 2 * L; i += 2)
                {
                    //My initial guess for nonlinear solver, a bit faster
                    
                    us_init_X[i]     = 1.5 * nodes[i / 2].velocity * v_sign[i / 2] - 0.5*us_init_X[i];
                    us_init_X[i + 1] = 1.5 * nodes[i / 2].lumen_sq - 0.5*us_init_X[i + 1];
                    
                    //Common initial guess for nonlinear solver
                  //  us_init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                  //  us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                    double wave_speed = 4 * (Math.Sqrt(Math.Sqrt(nodes[i / 2].lumen_sq)) * Math.Sqrt(beta_1[i / 2] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[i / 2]);
                    chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + wave_speed;
                    chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - wave_speed;
                }

                nl_system.solve(us_init_X, 1e-7, us_solution);

                double av_pressure = 0;
                double av_flux_in = 0;
                double av_flux_out = 0;
                double av_lumen_in = 0;
                double av_lumen_out = 0;

                for (int i = 0; i < 2 * L; i += 2)
                {                   

                    nodes[i / 2].velocity = us_solution[i] * v_sign[i / 2];
                    nodes[i / 2].lumen_sq = us_solution[i + 1];
                    nodes[i / 2].pressure = beta_1[i / 2] * (Math.Sqrt(nodes[i / 2].lumen_sq) - Math.Sqrt(nodes[i / 2].lumen_sq_0)) + GlobalDefs.DIASTOLIC_PRESSURE;

                    av_pressure += nodes[i / 2].pressure;
                    if (us_solution[i] >= 0)
                    {
                        av_flux_in += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_in += us_solution[i + 1];
                    }
                    else
                    {
                        av_flux_out += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_out += us_solution[i + 1];
                    }
                }

                core_node.pressure = av_pressure / L;
                core_node.lumen_sq = av_lumen_in;
                core_node.velocity = av_flux_in / av_lumen_in;
            }

            core_node.velocity = core_node.neighbours.Last().velocity;
            core_node.lumen_sq = core_node.neighbours.Last().lumen_sq;
#endif

        }

        protected double calcLumen_sq(int i, double pressure)
        {
            return (double)Math.Pow((pressure - GlobalDefs.DIASTOLIC_PRESSURE) / beta_1[i] + Math.Sqrt(lumen_sq_0[i]), 2);
        }

        protected MDFunction[] chrt_func;
        protected MDFunction mass_conservation_func;
        protected MDFunction[] energy_conservation_func;

        protected MDFunction[] funcs;

        protected double[] next_neighbours_pressure;
        protected VascularNode[] next_neighbours;
        protected int[] next_neighbours_v_sign;

        protected double[] prev_velocity;
        protected double[] wall_thickhess;
        protected double dt;

        protected double[] beta_1;
        protected double[] c_dst;
        protected double[] dX;
        protected double[] g_energy;

        protected double[] chrt_b;
        protected double[] chrt_f;

        protected bool[,] dep_matrix;

#if FAST
        protected NewtonSolver nl_system;
#endif

        unsafe protected double* us_init_X; 
        unsafe protected double* us_solution;
    }

    public class
        ViscoElasticKnot : Knot
    {
        public ViscoElasticKnot(Knot _knot, GetBetaFunction getElasticBeta)
            : base(_knot.core_node, _knot.current_time)
        {
            int L = nodes.GetLength(0);
            lumen_sq_old = new double[L];
            chrt_func = new MDFunction[L];
            energy_conservation_func = new MDFunction[L - 1];
            funcs = new MDFunction[2 * L];

            wall_thickhess = new double[L];
            lumen_sq_0 = new double[L];

            beta_1 = new double[L];
            chrt_b = new double[L];
            chrt_f = new double[L];
            c_dst = new double[L];
            dep_matrix = new bool[2 * L, 2 * L];
            prev_velocity = new double[L];
#if FAST
            nl_system = new NewtonSolver(2 * L);
#endif

            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    dep_matrix[i, j] = false;

            for (int i = 0; i < L; i++)
            {
                double R0 = Math.Sqrt(nodes[i].lumen_sq_0 / Math.PI);
                beta_1[i] = getElasticBeta(R0, nodes[i].elasticity) / nodes[i].lumen_sq_0;
                wall_thickhess[i] = GlobalDefs.getBoileauWallThickness(R0, nodes[i].elasticity);
                lumen_sq_0[i] = nodes[i].lumen_sq_0;
                prev_velocity[i] = nodes[i].velocity;
                chrt_b[i] = 0;
                chrt_f[i] = 0;               
                c_dst[i] = Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f /  GlobalDefs.BLOOD_DENSITY);
                lumen_sq_old[i] = lumen_sq_0[i];
            }

            for (int i = 0; i < L; i++)            
                pressure[i] = GlobalDefs.DIASTOLIC_PRESSURE;

            

            int count = 0;
            unsafe
            {
                for (int i = 0; i < L; i++)
                {
                    int I = i;
#if FAST
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v > 0)
                            return Math.Abs(v) + 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_b[I];
                    };
                    baseMDFunction f1 = new delegateMDFunc(f1_del);
                    nl_system.addFunc(f1);
                    nl_system.setDetMatrixEl(count, 2 * I, true);
                    nl_system.setDetMatrixEl(count, 2 * I + 1, true);
#endif
                    chrt_func[i] = delegate(double[] args) //v1,l1; v2,l2 ...
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v > 0)
                            return Math.Abs(v) + 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[I]) - chrt_b[I];
                    };

                    
                    funcs[count] = chrt_func[i];
                    dep_matrix[count, 2 * I] = true;
                    dep_matrix[count, 2 * I + 1] = true;

                    count++;
                }
            }

            unsafe
            {
#if FAST
                MDFunction_del f1_del = delegate(double* args)
                {
                    double summ_flux = 0;
                    for (int i = 0; i < L; i++)
                        summ_flux += args[0 + i * 2] * args[1 + i * 2];

                    return summ_flux;
                };
                baseMDFunction f1 = new delegateMDFunc(f1_del);
                nl_system.addFunc(f1);

                for (int i = 0; i < 2 * L; i++)
                    nl_system.setDetMatrixEl(count, i, true);                
#endif

                mass_conservation_func = delegate(double[] args)
                {
                    double summ_flux = 0;
                    for (int i = 0; i < L; i++)
                    {
                        double v = args[0 + i * 2];
                        double l = args[1 + i * 2];

                        summ_flux += v * l;
                    }
                    return summ_flux;
                };


                funcs[count] = mass_conservation_func;
                for (int i = 0; i < 2 * L; i++)                
                    dep_matrix[count, i] = true;
            };

            count++;

            unsafe
            {
                for (int i = 1; i < L; i++)
                {
                    int I = i;
#if FAST
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v0 = args[0];
                        double p0 = calcPressureV(0, args[1]);

                        double v = args[0 + I * 2];
                        double p = calcPressureV(I, args[1 + 2 * I]); 

                        return  GlobalDefs.BLOOD_DENSITY * (v0 * v0 - v * v) / 2 + p0 - p;
                    };
                    baseMDFunction f1 = new delegateMDFunc(f1_del);
                    nl_system.addFunc(f1);
                    nl_system.setDetMatrixEl(count, 0, true);
                    nl_system.setDetMatrixEl(count, 1, true);
                    nl_system.setDetMatrixEl(count, 2 * I, true);
                    nl_system.setDetMatrixEl(count, 2 * I + 1, true);
                   
#endif

                    energy_conservation_func[i - 1] = delegate(double[] args)
                    {
                        double v0 = args[0];
                        double p0 = calcPressureV(0, args[1]);

                        double v = args[0 + I * 2];
                        double p = calcPressureV(I, args[1 + 2 * I]);
                            
                        return  GlobalDefs.BLOOD_DENSITY * (v0 * v0 - v * v) / 2 + p0 - p;
                    };

                    
                    funcs[count] = energy_conservation_func[I - 1];

                    dep_matrix[count, 0] = true;
                    dep_matrix[count, 1] = true;
                    dep_matrix[count, 2 * I] = true;
                    dep_matrix[count, 2 * I + 1] = true;

                    count++;
                }

#if FAST
                us_init_X = (double*)Marshal.AllocHGlobal(2 * L * sizeof(double));
                us_solution = (double*)Marshal.AllocHGlobal(2 * L * sizeof(double));
                for (int i = 0; i < 2 * L; i += 2)
                {
                    us_init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                    us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                    nl_system.setDxVectorEl(i, 1e-12f);
                    nl_system.setDxVectorEl(i + 1, 1e-12f);
                }
#endif
            }

            dX = new double[2 * L];
            for (int i = 0; i < 2 * L; i += 2)
            {
                dX[i] = 1e-12f;
                dX[i + 1] = 1e-12f;
            }
        }
        

        unsafe public override void doCoupling(double dt)
        {
            curr_dt = dt;
            previous_time = current_time;
            current_time = current_time + dt;
            int L = nodes.GetLength(0);

            //if (core_node.id == 0)
                //L = nodes.GetLength(0);
#if SAFE
            double[] solution = new double[2 * L];
            double[] init_X = new double[2 * L];
            for (int i = 0; i < 2 * L; i += 2)
            {
                init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                init_X[i + 1] = nodes[i / 2].lumen_sq;

                double wave_speed = 4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[i / 2]);
                // Comment of nex two rows gives new algo for knot couplong
                chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + wave_speed;//4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[i/2]);
                chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - wave_speed;//4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f /  GlobalDefs.BLOOD_DENSITY) - c_dst[i/2]); 
            }

            solution = myMath.MDNewtonSolver(funcs, dep_matrix, init_X, 1e-6, dX);

            double av_pressure = 0;
            double av_flux_in = 0;
            double av_flux_out = 0;
            double av_lumen_in = 0;
            double av_lumen_out = 0;
            for (int i = 0; i < 2 * L; i += 2)
            {
                nodes[i / 2].velocity = solution[i] * v_sign[i / 2];
                nodes[i / 2].lumen_sq = solution[i + 1];
                nodes[i / 2].pressure = calcPressureV(i / 2, nodes[i / 2].lumen_sq);

                av_pressure += nodes[i / 2].pressure;
                if (solution[i] >= 0)
                {
                    av_flux_in += Math.Abs(solution[i + 1] * solution[i]);
                    av_lumen_in += solution[i + 1];
                }
                else
                {
                    av_flux_out += Math.Abs(solution[i + 1] * solution[i]);
                    av_lumen_out += solution[i + 1];
                }
            }
            core_node.velocity = core_node.neighbours.Last().velocity;
            core_node.lumen_sq = core_node.neighbours.Last().lumen_sq;
#endif

#if FAST1                        
            unsafe
            {
                for (int i = 0; i < 2 * L; i += 2)
                {
                    us_init_X[i    ] = 1.5*nodes[i / 2].velocity * v_sign[i / 2] - 0.5*us_init_X[i];
                    us_init_X[i + 1] = 1.5*nodes[i / 2].lumen_sq - 0.5*us_init_X[i + 1];
                    
                    double wave_speed = 4 * (Math.Sqrt(Math.Sqrt(nodes[i / 2].lumen_sq)) * Math.Sqrt(beta_1[i / 2] / 2.0f / GlobalDefs.BLOOD_DENSITY) - c_dst[i / 2]);
                    chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + wave_speed;
                    chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - wave_speed;
                }


                nl_system.solve(us_init_X, 1e-6, us_solution);


                double av_pressure = 0;
                double av_flux_in = 0;
                double av_flux_out = 0;
                double av_lumen_in = 0;
                double av_lumen_out = 0;

                for (int i = 0; i < 2 * L; i += 2)
                {
                    nodes[i / 2].velocity = us_solution[i] * v_sign[i / 2];
                    lumen_sq_old[i / 2] = nodes[i / 2].lumen_sq;
                    nodes[i / 2].lumen_sq = us_solution[i + 1];
                    nodes[i / 2].pressure = calcPressureV(i / 2, nodes[i / 2].lumen_sq);

                    av_pressure += nodes[i / 2].pressure;
                    if (us_solution[i] >= 0)
                    {
                        av_flux_in += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_in += us_solution[i + 1];
                    }
                    else
                    {
                        av_flux_out += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_out += us_solution[i + 1];
                    }
                }

                core_node.pressure = av_pressure / L;
                core_node.lumen_sq = av_lumen_in;
                core_node.velocity = av_flux_in / av_lumen_in;
            }


            core_node.velocity = core_node.neighbours.Last().velocity;
            core_node.lumen_sq = core_node.neighbours.Last().lumen_sq;
#endif
        }

        public override void reset()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                velocity[i] = 0;
                pressure[i] = GlobalDefs.DIASTOLIC_PRESSURE;                
                lumen_sq[i] = nodes[i].lumen_sq_0;
                lumen_sq_old[i] = nodes[i].lumen_sq_0;
                nodes[i].velocity = velocity[i];
                nodes[i].pressure = pressure[i];
                nodes[i].lumen_sq = nodes[i].lumen_sq_0;
            }
#if FAST
            unsafe
            {
                for (int i = 0; i < 2 * L; i += 2)
                {
                    us_init_X[i] = nodes[i / 2].velocity * v_sign[i / 2];
                    us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                }
            }
#endif
        }

        protected double calcPressureV(int i, double _lumen_sq)
        {
            double Gamma = 2.0 / 3.0 * Math.Sqrt(Math.PI) * wall_thickhess[i] * GlobalDefs.phi;
            return beta_1[i] * (Math.Sqrt(_lumen_sq) - Math.Sqrt(lumen_sq_0[i])) + GlobalDefs.DIASTOLIC_PRESSURE + Gamma / lumen_sq_0[i] / Math.Sqrt(_lumen_sq) * (_lumen_sq - lumen_sq_old[i]) / curr_dt;
        }

        protected double[] lumen_sq_old;
        double curr_dt;

        protected MDFunction[] chrt_func;
        protected MDFunction mass_conservation_func;
        protected MDFunction[] energy_conservation_func;

        protected MDFunction[] funcs;

        protected double[] next_neighbours_pressure;
        protected VascularNode[] next_neighbours;
        protected int[] next_neighbours_v_sign;

        protected double[] prev_velocity;
        protected double[] wall_thickhess;
        protected double dt;

        protected double[] beta_1;
        protected double[] c_dst;
        protected double[] dX;

        protected double[] chrt_b;
        protected double[] chrt_f;

        protected bool[,] dep_matrix;

#if FAST
        protected NewtonSolver nl_system;
        unsafe protected double* us_init_X;
        unsafe protected double* us_solution; 
#endif

       
    }


}

