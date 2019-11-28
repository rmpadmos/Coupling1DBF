#define FAST

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Runtime.InteropServices;

namespace TestModel
{
   

    public struct Vector3
    {
        public Vector3(double _x, double _y, double _z)
        { x = _x; y = _y; z = _z; }
        public double x, y, z;

        public static Vector3 operator *(Vector3 v1, double d)
        {
            return new Vector3(v1.x * d, v1.y * d, v1.z * d);
        }

        public static Vector3 operator -(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        }

        public static Vector3 operator -(Vector3 v1)
        {
            return new Vector3(-v1.x, -v1.y, -v1.z);
        }

        public static Vector3 operator +(Vector3 v1, Vector3 v2)
        {
            return new Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
        }

        static public double Distance(Vector3 v1, Vector3 v2)
        {
            return (double)Math.Sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z));
        }

        static public double Dot(Vector3 v1, Vector3 v2)
        {
            return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
        }

        static public double Length(Vector3 v1)
        {
            return Math.Sqrt((v1.x * v1.x + v1.y * v1.y + v1.z * v1.z));
        }

        public double Normilize()
        {
            double l = (double)Math.Sqrt(x * x + y * y + z * z);
            if (l == 0)
            {
                x = 0;
                y = 0;
                z = 0;
                return 0;
            }
            x = x / l;
            y = y / l;
            z = z / l;
            return l;
        }

    };

    public delegate double SimpleFunction (double   x );
    public delegate double GetBetaFunction(double   R0);
    public delegate double     MDFunction (double[] x );


    

    public class Thread
    {
        public Thread(List<Node> protothread)
        {
            

            int L = protothread.Count;
            nodes = (Node[])protothread.ToArray().Clone();
            v_sign = new int[L];

            velocity = new double[L];
            pressure = new double[L];
            lumen_sq = new double[L];
            flux     = new double[L];
            chrt     = new double[L];
       
            DefineSigns();
            current_time = 0;
        }

        unsafe public void nodes2state()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                velocity[i] = nodes[i].velocity * v_sign[i];
                pressure[i] = nodes[i].pressure;
                lumen_sq[i] = nodes[i].lumen_sq;
                chrt    [i] = nodes[i].chrt    ;
            }
        }

        public virtual void updateState()
        { }

        unsafe public void state2nodes()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                nodes[i].velocity = velocity[i] * v_sign[i];
                nodes[i].pressure = pressure[i];
                nodes[i].lumen_sq = lumen_sq[i];
                if (nodes[i].neighbours.Count < 3)
                    nodes[i].chrt = chrt[i];
                else
                    nodes[i].chrt = 0;
            }
        }

        public int getLength()
        {
            return nodes.GetLength(0);
        }

        public virtual void calcThread(double dt)
        {    
        }

        protected void DefineSigns()
        {
            int L = nodes.GetLength(0);
            Vector3 [] dir_vector1 = new Vector3[L];
            Vector3 [] dir_vector2 = new Vector3[L];
                            v_sign = new int    [L];
            for (int i = 0; i < L - 1; i++)
            {
                dir_vector1[i] = (nodes[i + 1].position - nodes[i].position);
            }
            for (int i = 0; i < L; i++)
            {
                nodes[i].defDirVector();
                dir_vector2[i] = nodes[i].dir_vector;
            }

            if(L>1)
                dir_vector1[L - 1] = nodes[L - 1].position - nodes[L - 2].position;
            else
                dir_vector1[L - 1] = dir_vector2[L - 1];
                        
            for (int i = 0; i < L; i++)
                v_sign[i] = Math.Sign(Vector3.Dot(dir_vector1[i], dir_vector2[i]));
        }

        public void WriteThreadToFile(string file_path)
        {
            int id_count = 0;
            string outText = "Name:";
            outText += "System_0" + "\n";
            outText += "Coordinates:\n";
            foreach (var n in nodes)
            {
                n.id = id_count;
                outText += n.id + " X:" + n.position.x.ToString("F8") + " Y:" + n.position.y.ToString("F8") + " Z:" + n.position.z.ToString("F8") + " R:" + Math.Sqrt(n.lumen_sq_0 / Math.PI).ToString("F8") + "\n";
                id_count++;
            }
            outText += "\nBonds:\n";
            foreach (var n in nodes)
            {
                outText += n.id + " ";
                foreach (var nn in n.neighbours)
                    if (nodes.Contains(nn))
                    {
                        outText += nn.id + " ";
                    }
                outText += "\n";
            }
            System.IO.File.WriteAllText(file_path, outText);
        }

        public virtual void setClot(Node nd, float degree)
        { }

        public virtual void removeClot(Node nd)
        { }


        public Node[] nodes;
        public int[] v_sign;

        public double[] velocity;
        public double[] pressure;
        public double[] lumen_sq;
        public double[] flux;
        public double[] chrt;

        public double current_time;
    }


    public class ElasticThread : Thread 
    {

        public ElasticThread(Thread _protothread, GetBetaFunction getElsticBeta):base(_protothread.nodes.ToList())
        {            
            int L = nodes.GetLength(0);
            beta_1      = new double[L];
            lumen_sq_0  = new double[L];

            for (int i = 0; i < L; i++)
            {
                double R0 = Math.Sqrt(nodes[i].lumen_sq_0 / Math.PI);                
                beta = getElsticBeta(R0);
                beta_1    [i] = beta / nodes[i].lumen_sq_0;
                lumen_sq_0[i] = nodes[i].lumen_sq_0;
                lumen_sq  [i] = nodes[i].lumen_sq_0;
                pressure  [i] = Node.diastolic_pressure;                                
            }

            unsafe
            {
                double* dZ = (double*)Marshal.AllocHGlobal((L-1) * sizeof(double));
                for (int i = 0; i < L-1; i++)
                    dZ[i] = Vector3.Distance(nodes[i + 1].position, nodes[i].position);

                thread_solver = new McCormackThread(L, dZ, Node.density, Node.viscosity); 

                for (int i = 0; i < L; i++)
                {
                    double bt = beta_1[i];
                    double lm = lumen_sq_0[i];
                    _1DFunction_del pressure_func_del = delegate(double x)
                    {
                        return Node.diastolic_pressure + (bt * (Math.Sqrt(x) - Math.Sqrt(lm)));
                    };

                    delegate1DFunc f = new delegate1DFunc(pressure_func_del);

                    thread_solver.addFunc(f, i);                    
                }                
            }

            clotes = new List<Clot>();
        }

        public override void setClot(Node nd, float degree)
        {
            int nd_id = 0;
            for (int i = 0; i < nodes.GetLength(0); i++)
                if (nodes[i] == nd)
                { nd_id = i; break; }

            if (clotes.FindIndex(x => x.c_nd == nd) == -1)
            {

                Clot clt = new Clot();
                clt.c_nd = nd;

                clt.normal_lumen = new List<double>();
                clt.nd = new List<Node>();

                int L = 2;

                for(int i=-L; i<L; i++)
                    if(nd_id+i>0&&nd_id+i<nodes.GetLength(0))
                    {
                        clt.nd.Add(nodes[nd_id + i]);
                        clt.normal_lumen.Add(lumen_sq_0[nd_id + i]);
                    }

                clt.tgr_degree = degree;
                clt.curr_degree = 0;
                clt.thread_id = nd_id;
                clotes.Add(clt);
            }
        }

        public override void removeClot(Node nd)
        {
            int clt_id = clotes.FindIndex(x => x.c_nd == nd);
            if (clt_id != -1)
            {   
                Clot clt = clotes[clt_id];
                clt.tgr_degree = 0.0;
                clotes[clt_id] = clt;
            }
        }

        public override void updateState()
        {
            for (int i = 0; i < clotes.Count; i++)
            {
                Clot clt = clotes[i];
                int id = clt.thread_id;
                if (clt.tgr_degree - clt.curr_degree > 0.01)
                {
                    clt.curr_degree += 0.01f;

                    int L = 2;
                    int nd_id = clt.thread_id;
                    for (int j = -L; j < L; j++)
                        if (nd_id + j > 0 && nd_id + j < nodes.GetLength(0))
                        {
                            double degree = clt.curr_degree *(L-Math.Abs(j)+1)/(L+1);
                            clt.nd[j + L].lumen_sq_0 = clt.normal_lumen[j + L] * (1 - degree);
                            lumen_sq_0[id + j] = clt.nd[j + L].lumen_sq_0;
                            beta_1[id + j] = this.beta/clt.nd[j + L].lumen_sq_0;
                        }
                    clotes[i] = clt;
                }
                else
                {
                    if (clt.tgr_degree == 0.0)
                    {
                        clotes.RemoveAt(i);
                        break;
                    }
                }
            }
        }

        unsafe public override void calcThread(double dt)
        {
            current_time += dt;
            double dz = 0;
           
           int L = nodes.GetLength(0);

           this.nodes2state();

          


       //    double[] velocity = (double[])velocity.Clone();
       //    double[] lumen_sq = (double[])lumen_sq.Clone();
       //    double[] pressure = (double[])pressure.Clone();

       /*    try
           {
               this.nodes.First(x => x.id == 2801);
               this.nodes[0].id = this.nodes[0].id;
           }
           catch { }*/


           unsafe 
           {
               fixed (double* v_ptr = &velocity[0])
               {
                   fixed (double* lum_ptr = &lumen_sq[0])
                   {
                       fixed (double* p_ptr = &pressure[0])
                       {
                           thread_solver.calc(v_ptr, lum_ptr, p_ptr, dt);
                       }
                   }
               }
           }

            /*
           double[] velocity_pred = (double[])velocity.Clone();
           double[] lumen_sq_pred = (double[])lumen_sq.Clone();
           double[] pressure_pred = (double[])pressure.Clone();

           for (int i = 1; i < L - 1; i++)
           {
               dz = Vector3.Distance(nodes[i].position, nodes[i - 1].position);
               velocity_pred[i] = velocity[i] - dt / dz * ((velocity[i] * velocity[i] / 2 + pressure[i] / Node.density) - (velocity[i - 1] * velocity[i - 1] / 2 + pressure[i - 1] / Node.density));
               lumen_sq_pred[i] = lumen_sq[i] - dt / dz * (lumen_sq[i] * velocity[i] - lumen_sq[i-1] * velocity[i-1]);
               pressure_pred[i] = calcPressure(lumen_sq_pred[i],i);
           }           

           for (int i = 1; i < L - 1; i++)
           {
               dz = Vector3.Distance(nodes[i].position, nodes[i + 1].position);
               velocity[i] = (velocity[i] + velocity_pred[i]) / 2 - dt / dz / 2 * ((velocity_pred[i + 1] * velocity_pred[i + 1] / 2 + pressure_pred[i + 1] / Node.density) - (velocity_pred[i] * velocity_pred[i]/ 2 + pressure_pred[i] / Node.density));
               velocity   [i] =  velocity[i] - 1.0 / Node.density / lumen_sq[i] * (22 * Node.viscosity * Math.PI * velocity[i]) * dt;
               lumen_sq   [i] = (lumen_sq[i] + lumen_sq_pred[i]) / 2 - dt / dz / 2 * (lumen_sq_pred[i + 1] * velocity_pred[i + 1] - lumen_sq_pred[i] * velocity_pred[i]);
               pressure   [i] = calcPressure(lumen_sq[i],   i);               
           }

           if (L > 1)
           {
               velocity[0]     = velocity[1];
               velocity[L - 1] = velocity[L - 2];

               pressure[0] = pressure[1];
               pressure[L - 1] = pressure[L - 2];

               lumen_sq[0] = lumen_sq[1];
               lumen_sq[L - 1] = lumen_sq[L - 2];
           }            
            */
           this.state2nodes();
        }

        protected double calcPressure(double _lumen_sq, int i)
        {
            return Node.diastolic_pressure + (beta_1[i] * (Math.Sqrt(_lumen_sq) - Math.Sqrt(lumen_sq_0[i])));
        }


        protected double calcPressure(int i)
        {
            return Node.diastolic_pressure + (beta_1[i] * (Math.Sqrt(lumen_sq[i]) - Math.Sqrt(lumen_sq_0[i])));
        }


        protected double calcLumen_sq(int i)
        {   
            return (double)Math.Pow((pressure[i] - Node.diastolic_pressure) / beta_1[i] + Math.Sqrt(lumen_sq_0[i]), 2);
        }
         
        protected double    young_modulus   ;
        protected double    wall_thickhess  ;
        protected double    beta            ;
        protected double[]  beta_1          ;
        protected double[]  lumen_sq_0      ;

        McCormackThread     thread_solver;

        const double h_a = 0.2802;
        const double h_b = -505.3; //m^-1
        const double h_c = 0.1324;
        const double h_d = -11.14; //m^-1 

        protected List<Clot> clotes;
    }

    public class Knot
    {        
        public Knot(Node _core_node, double start_time)
        {
            core_node = _core_node;
            nodes = (Node[])core_node.neighbours.ToArray().Clone();

            int L = nodes.GetLength(0);

              v_sign = new    int[L];
            velocity = new double[L];
            pressure = new double[L];
            lumen_sq = new double[L];

            for (int i = 0; i < core_node.neighbours.Count; i++ )
                lumen_sq[i] = core_node.neighbours[i].lumen_sq_0;

            DefineSigns();

            current_time  = start_time;
            previous_time = current_time - 1e-3;
        }

        public void DefineSigns()
        {
            int L = nodes.GetLength(0);
            Vector3[] dir_vector1 = new Vector3[L];
            Vector3[] dir_vector2 = new Vector3[L];
            v_sign = new int[L];
            for (int i = 0; i < L; i++)
            {
                dir_vector1[i] =  core_node.position - nodes[i].position;
                dir_vector2[i] =  nodes[i].dir_vector;
            }
            
            for (int i = 0; i < L; i++)
                v_sign[i] = Math.Sign(Vector3.Dot(dir_vector1[i], dir_vector2[i]));
        }

        public virtual void doCoupling(double dt){}

        public virtual void holdChtr(double dt) { }

        public Node[]       nodes;
        public Node     core_node;


        public int   []   v_sign;
        public double[] velocity;
        public double[] pressure;
        public double[] lumen_sq;

        public  double  current_time;
        public  double previous_time;
    }

    public class 
        StandartKnot : Knot
    {
        public StandartKnot(Knot _knot, GetBetaFunction getElasticBeta): base(_knot.core_node, _knot.current_time)
        {
            int L = nodes.GetLength(0);
            chrt_func                = new MDFunction[L];
            energy_conservation_func = new MDFunction[L - 1];
            funcs                    = new MDFunction[2 * L];
            

            beta_1 = new double[L];
            chrt_b = new double[L];
            chrt_f = new double[L];
            c_dst  = new double[L];
            dep_matrix = new bool[2*L, 2*L];
            prev_velocity = new double[L];
       
            nl_system = new NewtonSolver(2*L);


            for (int i = 0; i < L; i++)
                for (int j = 0; j < L; j++)
                    dep_matrix[i, j] = false;

            for (int i = 0; i < L; i++)
            {
                double R0 = Math.Sqrt(nodes[i].lumen_sq_0 / Math.PI);
                beta_1[i] = getElasticBeta(R0) / nodes[i].lumen_sq_0;
                prev_velocity[i] = nodes[i].velocity;             
                chrt_b[i] = 0;//-(4 * Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f / Node.density));
                chrt_f[i] = 0;// (4 * Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f / Node.density));               
                 c_dst[i] = Math.Pow(nodes[i].lumen_sq_0, 0.25f) * Math.Sqrt(beta_1[i] / 2.0f / Node.density);
            }   

            int count = 0;
            unsafe{
                for (int i = 0; i < L; i++)
                {
                    int I = i;
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v > 0)
                            return Math.Abs(v) + 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f / Node.density) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Sqrt(Math.Sqrt(l)) * Math.Sqrt(beta_1[I] / 2.0f / Node.density) - c_dst[I]) - chrt_b[I]; 
                    };
                    baseMDFunction f1 = new delegateMDFunc(f1_del);

                    chrt_func[i] = delegate(double[] args) //v1,l1; v2,l2 ...
                    {
                        double v = args[0 + I * 2];
                        double l = args[1 + I * 2];

                        if (v>0)
                            return Math.Abs(v) + 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f / Node.density) - c_dst[I]) - chrt_f[I];
                        else
                            return Math.Abs(v) - 4 * (Math.Pow(l, 0.25f) * Math.Sqrt(beta_1[I] / 2.0f / Node.density) - c_dst[I]) - chrt_b[I];                    
                    };

                    nl_system.addFunc(f1);
                    funcs[count] = chrt_func[i];

                    dep_matrix[count, 2 * I]     = true;
                    dep_matrix[count, 2 * I + 1] = true;
                    nl_system.setDetMatrixEl(count, 2 * I  ,true);
                    nl_system.setDetMatrixEl(count, 2 * I+1,true);

                    count++;                
                }
            }

            unsafe
            {
                MDFunction_del f1_del = delegate(double* args)
                {
                    double summ_flux = 0;
                    for (int i = 0; i < L; i++)     
                        summ_flux += args[0 + i * 2] * args[1 + i * 2];
                 
                    return summ_flux;
                };
                baseMDFunction f1 = new delegateMDFunc(f1_del);


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
                {
                    dep_matrix[count, i] = true;
                    nl_system.setDetMatrixEl(count, i, true);
                }

                nl_system.addFunc(f1);
            };

            count++;

            unsafe
            {
                for (int i = 1; i < L; i++)
                {
                    int I = i;
                    MDFunction_del f1_del = delegate(double* args)
                    {
                        double v0 = args[0];
                        double p0 = beta_1[0] * (Math.Sqrt(args[1]) - Math.Sqrt(nodes[0].lumen_sq_0)) + Node.diastolic_pressure;

                        double v = args[0 + I * 2];
                        double p = beta_1[I] * (Math.Sqrt(args[1 + 2 * I]) - Math.Sqrt(nodes[I].lumen_sq_0)) + Node.diastolic_pressure;
                        return Node.density * (v0 * v0 - v * v) / 2 + p0 - p;
                    };
                    baseMDFunction f1 = new delegateMDFunc(f1_del);

                    energy_conservation_func[i - 1] = delegate(double[] args)
                    {
                        double v0 = args[0];
                        double p0 = beta_1[0] * (Math.Sqrt(args[1]) - Math.Sqrt(nodes[0].lumen_sq_0)) + Node.diastolic_pressure;

                        double v = args[0 + I * 2];
                        double p = beta_1[I] * (Math.Sqrt(args[1 + 2 * I]) - Math.Sqrt(nodes[I].lumen_sq_0)) + Node.diastolic_pressure;
                        return Node.density * (v0 * v0 - v * v) / 2 + p0 - p;
                    };

                    nl_system.addFunc(f1);
                    funcs[count] = energy_conservation_func[I - 1];

                    dep_matrix[count, 0] = true;
                    dep_matrix[count, 1] = true;
                    nl_system.setDetMatrixEl(count, 0, true);
                    nl_system.setDetMatrixEl(count, 1, true);

                    dep_matrix[count, 2 * I] = true;
                    dep_matrix[count, 2 * I + 1] = true;

                    nl_system.setDetMatrixEl(count, 2 * I    , true);
                    nl_system.setDetMatrixEl(count, 2 * I + 1, true);

                    count++;
                }
            }

            dX = new double[2 * L];
            for (int i = 0; i < 2 * L; i += 2)
            {
                dX[i]     = 1e-12f;
                dX[i + 1] = 1e-12f;
                nl_system.setDxVectorEl(i    , 1e-12f);
                nl_system.setDxVectorEl(i + 1, 1e-12f);
            }
        }

        public override void holdChtr(double dt)
        {
            int L = nodes.GetLength(0);
            for (int i=0; i<L; i++)
            {
                if (nodes[i].velocity * v_sign[i] > 0)
                {
                    Node node_1 = nodes[i].neighbours.Find(x => x != core_node);
                    double c0 = Math.Sqrt(beta_1[i] / 2 / Node.density) * Math.Pow(nodes[i].lumen_sq, 0.25);                    
                    double c1 = Math.Sqrt(beta_1[i] / 2 / Node.density) * Math.Pow(node_1.lumen_sq, 0.25);
                    double lambda_f = c0 + Math.Abs(nodes[i].velocity);
                    double chrt_f_0 = Math.Abs(nodes[i].velocity) + 4 * (c0 - c_dst[i]);
                    double chrt_f_1 = Math.Abs(node_1.velocity  ) + 4 * (c1 - c_dst[i]);                   
                    double fact = lambda_f * dt / Vector3.Length(nodes[i].position - node_1.position);
                    chrt_f[i] = (chrt_f_0 * (1 - fact) + chrt_f_1 * fact);
                }
                else
                {
                    Node node_1 = nodes[i].neighbours.Find(x => x != core_node);
                    double c0 = Math.Sqrt(beta_1[i] / 2 / Node.density) * Math.Pow(nodes[i].lumen_sq, 0.25);
                    double c1 = Math.Sqrt(beta_1[i] / 2 / Node.density) * Math.Pow(  node_1.lumen_sq, 0.25);
                    double lambda_b = c0 - Math.Abs(nodes[i].velocity);
                    double chrt_b_0 = Math.Abs(nodes[i].velocity) - 4 * (c0 - c_dst[i]);
                    double chrt_b_1 = Math.Abs(node_1.velocity  ) - 4 * (c1 - c_dst[i]);
                    double fact = lambda_b * dt / Vector3.Length(nodes[i].position - node_1.position);
                    chrt_b[i] = (chrt_b_0 * (1 - fact) + chrt_b_1 * fact);
                }
            }
        }

        unsafe public override void doCoupling(double dt)
        {        
            current_time  = current_time + dt;
            previous_time = current_time;
            
            int L = nodes.GetLength(0);
           
            
#if SAFE            
            double[] solution = new double[2 * L];
            double[] init_X   = new double[2 * L];     
            for(int i=0; i<2*L; i+=2)
            {
                init_X[i  ] = nodes[i/2].velocity*v_sign[i/2];
                init_X[i+1] = nodes[i/2].lumen_sq;
                double wave_speed = 4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f / Node.density) - c_dst[i/2]);
                // Comment of nex two rows gives new algo for knot couplong
                chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + 4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f / Node.density) - c_dst[i/2]);
                chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - 4 * (Math.Pow(nodes[i / 2].lumen_sq, 0.25f) * Math.Sqrt(beta_1[i / 2] / 2.0f / Node.density) - c_dst[i/2]); 
            }
            
            solution = myMath.MDNewtonSolver(funcs, dep_matrix, init_X, 1e-6, dX);

            for (int i = 0; i < 2 * L; i += 2)
            {
                nodes[i / 2].velocity = solution[i] * v_sign[i / 2];
                nodes[i / 2].lumen_sq = solution[i + 1];
                nodes[i / 2].pressure = beta_1[i / 2] * (Math.Sqrt(nodes[i / 2].lumen_sq) - Math.Sqrt(nodes[i / 2].lumen_sq_0)) + Node.diastolic_pressure;
            }
#endif

#if FAST
            unsafe
            {
                double* us_init_X   = (double*)Marshal.AllocHGlobal(2* L * sizeof(double));
                double* us_solution = (double*)Marshal.AllocHGlobal(2* L * sizeof(double));

                for (int i = 0; i < 2 * L; i += 2)
                {
                //    double wave_speed = 4 * Math.Sqrt(Math.Sqrt(nodes[i / 2].lumen_sq)) * Math.Sqrt(beta_1[i / 2] / 2.0f / Node.density) - c_dst[i / 2];

                    us_init_X[i    ] = nodes[i / 2].velocity * v_sign[i / 2];
                    us_init_X[i + 1] = nodes[i / 2].lumen_sq;
                    double wave_speed = 4 * (Math.Sqrt(Math.Sqrt(nodes[i / 2].lumen_sq)) * Math.Sqrt(beta_1[i / 2] / 2.0f / Node.density) - c_dst[i / 2]);
                    chrt_f[i / 2] = Math.Abs(nodes[i / 2].velocity) + wave_speed;
                    chrt_b[i / 2] = Math.Abs(nodes[i / 2].velocity) - wave_speed; 
                }

                nl_system.solve(us_init_X, 1e-6, us_solution);

               // for (int i = 0; i < 2 * L; i++)
               //     solution[i] = us_solution[i];
              /*  double av_pressure = 0;
                double av_flux_in = 0;
                double av_flux_out = 0;
                double av_lumen_in = 0;
                double av_lumen_out = 0;*/

                for (int i = 0; i < 2 * L; i += 2)
                {
                    nodes[i / 2].velocity = us_solution[i] * v_sign[i / 2];
                    nodes[i / 2].lumen_sq = us_solution[i + 1];
                    nodes[i / 2].pressure = beta_1[i / 2] * (Math.Sqrt(nodes[i / 2].lumen_sq) - Math.Sqrt(nodes[i / 2].lumen_sq_0)) + Node.diastolic_pressure;

               /*     av_pressure += nodes[i / 2].pressure;
                    if (us_solution[i] >= 0)
                    {
                        av_flux_in += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_in += us_solution[i + 1];
                    }
                    else
                    {
                        av_flux_out += Math.Abs(us_solution[i + 1] * us_solution[i]);
                        av_lumen_out += us_solution[i + 1];
                    }*/
                }

                Marshal.FreeHGlobal((IntPtr)us_solution);
                Marshal.FreeHGlobal((IntPtr)us_init_X);

            }


            core_node.velocity = core_node.neighbours.Last().velocity;
            core_node.lumen_sq = core_node.neighbours.Last().lumen_sq;    

            /*


            double  av_pressure = 0;
            double  av_flux_in = 0;
            double  av_flux_out = 0;
            double  av_lumen_in = 0;
            double  av_lumen_out = 0;
            for (int i = 0; i < 2 * L; i += 2)
            {
                nodes[i / 2].velocity = solution[i  ] * v_sign[i/2];
                nodes[i / 2].lumen_sq = solution[i+1];
                nodes[i / 2].pressure = beta_1  [i/2] * (Math.Sqrt(nodes[i/2].lumen_sq) - Math.Sqrt(nodes[i/2].lumen_sq_0)) + Node.diastolic_pressure;              

                av_pressure += nodes[i / 2].pressure;
                if (solution[i] >= 0)
                {
                    av_flux_in += Math.Abs(solution[i + 1] * solution[i]);
                    av_lumen_in += solution[i + 1];
                }
                else
                {
                    av_flux_out  += Math.Abs(solution[i + 1] * solution[i]);
                    av_lumen_out += solution[i + 1];
                }
            }*/
#endif
        }

        protected MDFunction[] chrt_func;
        protected MDFunction   mass_conservation_func;
        protected MDFunction[] energy_conservation_func;

        protected MDFunction[] funcs;

        protected double[] next_neighbours_pressure;
        protected Node  [] next_neighbours;
        protected int   [] next_neighbours_v_sign;
        
        protected double[] prev_velocity;
        protected double   dt;

        protected double[] beta_1;
        protected double[] c_dst;
        protected double[] dX;

        protected double[] chrt_b;
        protected double[] chrt_f;

        protected   bool[,] dep_matrix;

        protected NewtonSolver nl_system;
    }



    public delegate int NodeFilter(Node node);
    public delegate double getFloatValueDelegate (Node node);
    public delegate void  setFloatValueDelegate (Node node, double value);

    public delegate void  setBoolValueDelegate  (Node node, bool value);
    public delegate bool  getBoolValueDelegate  (Node node);

    public class Node
    {

        public static double Wall_time
        {
            set { wall_time = value; }
            get { return wall_time; }
        }

        public static void ResetWallTime(double start_time, double _dt)
        {
            wall_time = start_time;
            dt = _dt;
        }

        public Node(int _id, Vector3 _position, double _rad)
        {
            id = _id;
            position = _position;
            neighbours = new List<Node>();
            lumen_sq_0 = (double)Math.PI * _rad * _rad;
            lumen_sq   = lumen_sq_0;
            pressure = Node.diastolic_pressure;
            velocity = 0;
            f_id2value = new Dictionary<int, double>();
            b_id2value = new Dictionary<int, bool>();
        }

        public virtual void defDirVector()
        {
            if (neighbours.Count < 3)
                dir_vector = neighbours.First().position - position;
            else
                dir_vector = new Vector3(0, 0, 0);
        }

        public int addNeighbours(Node[] _neighbours)
        {
            int L1 = neighbours.Count();
            neighbours.AddRange(_neighbours.ToList());
            neighbours = neighbours.Distinct().ToList<Node>();
            neighbours.RemoveAll(x => x == this);

            return neighbours.Count - L1;
        }

        public int excludeNeighbours(Node[] _neighbours)
        {
            int L1 = neighbours.Count;
            foreach (var n in _neighbours)
                excludeNeighbour(n);

            return L1 - neighbours.Count;
        }

        public void addNeighbour(Node neighbour)
        {
            neighbours.Remove(neighbour);
            neighbours.Add(neighbour);
        }

        public bool excludeNeighbour(Node neighbour)
        {
            return neighbours.Remove(neighbour);
        }

        public static void newFloatValueLayer(out getFloatValueDelegate get_del, out setFloatValueDelegate set_del)
        {
            Random r = new Random();

            int c = r.Next(int.MaxValue); ;
            while (f_dlg2id.ContainsValue(c))
                c = r.Next(int.MaxValue);

            get_del = delegate(Node node)
            {
                Node n = (Node)node;
                if (n.f_id2value.ContainsKey(c))
                    return n.f_id2value[c];
                else
                {
                    n.f_id2value.Add(c, double.MaxValue);
                    Node.id2nodelist[c].Add(node);
                }
                return double.MaxValue;
            };

            f_dlg2id.Add(get_del, c);

            set_del = delegate(Node node, double val)
            {
                Node n = (Node)node;
                if (n.f_id2value.ContainsKey(c))
                    n.f_id2value[c] = val;
                else
                {
                    n.f_id2value.Add(c, val);
                    Node.id2nodelist[c].Add(node);
                }
            };

            id2nodelist.Add(c, new List<Node>());
        }

        public static void newBoolValueLayer(out getBoolValueDelegate get_del, out setBoolValueDelegate set_del)
        {
            Random r = new Random();

            int c = r.Next(int.MaxValue); ;
            while (b_dlg2id.ContainsValue(c))
                c = r.Next(int.MaxValue);

            get_del = delegate(Node node)
            {
                Node n = (Node)node;
                if (n.b_id2value.ContainsKey(c))
                    return n.b_id2value[c];
                else
                {
                    n.b_id2value.Add(c, false);
                    Node.id2nodelist[c].Add(node);
                }
                return false;
            };

            b_dlg2id.Add(get_del, c);

            set_del = delegate(Node node, bool val)
            {
                Node n = (Node)node;
                if (n.b_id2value.ContainsKey(c))
                    n.b_id2value[c] = val;
                else
                {
                    n.b_id2value.Add(c, val);
                    Node.id2nodelist[c].Add(node);
                }
            };

            id2nodelist.Add(c, new List<Node>());
        }

        private Dictionary<int, double> f_id2value;
        private Dictionary<int, bool > b_id2value;
        private static Dictionary<getFloatValueDelegate, int> f_dlg2id = new Dictionary<getFloatValueDelegate, int>();
        private static Dictionary<getBoolValueDelegate, int> b_dlg2id = new Dictionary<getBoolValueDelegate, int>();
        private static Dictionary<int, List<Node>> id2nodelist = new Dictionary<int, List<Node>>();


        public static void terminateFloatValueLayer(ref getFloatValueDelegate get_del)
        {
            int c = f_dlg2id[get_del];
            foreach (Node n in id2nodelist[c])
                n.f_id2value.Remove(c);
            id2nodelist.Remove(c);
            f_dlg2id.Remove(get_del);
            get_del = null;
        }

        public static void terminateBoolValueLayer(ref getBoolValueDelegate get_del)
        {
            int c = b_dlg2id[get_del];
            foreach (Node n in id2nodelist[c])
                n.b_id2value.Remove(c);
            id2nodelist.Remove(c);
            b_dlg2id.Remove(get_del);
            get_del = null;
        }

        public double velocity { get; set; }
        public double lumen_sq { get; set; }
        public double pressure { get; set; }
        public double chrt     { get; set; }
        
        public static double diastolic_pressure = 10.9e+3f; // Pa
        public static double density            = 1060;     // kg/m3
        public static double viscosity          = 3.5e-3f;  // Pa*s
        public static double pulse_time_est     = 1.1f;     // s

        public int                id;
        public Vector3    dir_vector;        
        public List<Node> neighbours;
        public Vector3      position;
        public double     lumen_sq_0;
        
        protected static double wall_time;
           public static double        dt;
    };

    public class BoundaryCondition
    {
        public BoundaryCondition(Node _core_node, double start_time)
        {
            core_node       = _core_node;
            current_time    = start_time;
            previous_time   = current_time - 1e-3;
            v_sign          = new int[2];
            DefineSign();
        }

        public      virtual     void     doBC   (double dt)  
        {}
        public      virtual     double    getTime()           
        { return current_time; }


        protected virtual    void     DefineSign()     
        {
            Vector3 dir_vector1 = new Vector3();
            Vector3 dir_vector2 = new Vector3();
            

            dir_vector1 = core_node.neighbours.Last().position - core_node.position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));

            dir_vector2 = core_node.neighbours.Last().dir_vector;
            v_sign[1] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));
        }        

        public           double       current_time   ;
        public           double       previous_time  ;
        public             int[]       v_sign        ;
        public            Node       core_node       ;
    };

    public class InletPressure: BoundaryCondition
    {
        public InletPressure(BoundaryCondition BC, TableFunction _pressure, double _beta): base(BC.core_node, BC.current_time)
        {
            this.core_node      =    base.core_node    ;
            this.current_time   =    base.current_time ;
            this.previous_time  =    base.previous_time;
            this.v_sign         =    base.v_sign;
            this.previous_velocity = 0;

            pressure_on_time = _pressure;

            beta_1 = _beta / core_node.lumen_sq_0;

            if (core_node.neighbours.Find(x => x.neighbours.Count <= 2) != null)
            {
                core_node.lumen_sq_0 = core_node.neighbours.Find(x => x.neighbours.Count <= 2).lumen_sq_0;
                core_node.lumen_sq   = core_node.lumen_sq_0;
                core_node.chrt       = core_node.neighbours.Find(x => x.neighbours.Count <= 2).chrt;
            }

            core_node.chrt = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
            this.chrt = core_node.chrt;
        }

        public override void doBC(double dt)
        {
            if (true)
            {
                current_time = current_time + dt;
                previous_time = current_time;
                previous_velocity = current_velocity;

                double pressure = pressure_on_time(current_time);
                double inlet_lumen = Math.Pow((pressure - Node.diastolic_pressure) / beta_1 + Math.Sqrt(core_node.lumen_sq_0),2);
                double neighbour_chrt = core_node.neighbours.Last().velocity - 4 * Math.Pow(core_node.neighbours.Last().lumen_sq,0.25)*Math.Sqrt(beta_1/Node.density/2.0f);
                double U = neighbour_chrt + 4 * Math.Pow(inlet_lumen, 0.25) * Math.Sqrt(beta_1 / Node.density / 2.0f);        

                core_node.velocity = U * v_sign[0];
                core_node.pressure = pressure;
                core_node.lumen_sq = inlet_lumen;
                chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
                core_node.chrt = chrt;

                current_velocity = U;               
            }
        }

        protected TableFunction pressure_on_time;
        protected SimpleFunction flux_function;

        protected double                beta_1;
        protected double                Q_min ;
        protected double                Q_max ;

        protected double     previous_velocity;
        protected double      current_velocity;
        protected double                  chrt;

        protected double    diastolic_pressure;
        protected double     sistolic_pressure;

        protected Queue<double> pressure_hist;

        protected double         pulse_time_interval;
        protected double   flux_minmax_time_interval;

        private float sample_dt;
    };

    public class InletFlux : BoundaryCondition
    {
        public InletFlux(BoundaryCondition BC, TableFunction _flux, GetBetaFunction getElasticBeta)
            : base(BC.core_node, BC.current_time)
        {
            this.core_node = base.core_node;
            this.current_time = base.current_time;
            this.previous_time = base.previous_time;
            this.v_sign = base.v_sign;
            this.previous_velocity = 0;

            base_flux_on_time = _flux;
            base_period       =   1.0f;

            flux_on_time = _flux;

            

            if (core_node.neighbours.Find(x => x.neighbours.Count <= 2) != null)
            {
                core_node.lumen_sq_0 = core_node.neighbours.Find(x => x.neighbours.Count <= 2).lumen_sq_0;
                core_node.lumen_sq = core_node.lumen_sq_0;
                core_node.chrt = core_node.neighbours.Find(x => x.neighbours.Count <= 2).chrt;
            }

            double R0 = Math.Sqrt(core_node.lumen_sq_0 / Math.PI);
            beta_1 = getElasticBeta(R0) / core_node.lumen_sq_0;
            flux_function = delegate(double A)
            {
                return (A * core_node.chrt - 4 * Math.Pow(A, 1.25f) * Math.Sqrt(beta_1 / 2.0 / Node.density) - flux_on_time(current_time)) * 1000;
            };

            core_node.chrt = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
            this.chrt = core_node.chrt;


            sample_dt = 0.05f;
            double t_min = 0;
            double t_max = 0;
            Q_min = flux_on_time(current_time);
            Q_max = Q_min;
            for (int i = 1; i < Node.pulse_time_est / sample_dt; i++)
            {
                double val = flux_on_time(current_time + i * sample_dt);
                if (val > Q_max)
                {
                    Q_max = val;
                    t_max = i * sample_dt;
                }
                else
                    if (val < Q_min)
                    {
                        Q_min = val;
                        t_min = i * sample_dt;
                    }
            }
            flux_minmax_time_interval = Math.Abs(t_max - t_min);

            pulse_time_interval = 0;
            pressure_hist = new Queue<double>();
        }

        public override void doBC(double dt)
        {
            if (true)
            {
                current_time = current_time + dt;
                previous_time = current_time;
                previous_velocity = current_velocity;


                if (current_time > 23.0)
                    previous_time = current_time;

                double flux = flux_on_time(current_time);
                double inlet_lumen = core_node.neighbours.Last().lumen_sq;//myMath.NewtonSolver(flux_function, core_node.lumen_sq_0, 1e-10f, core_node.lumen_sq_0 * 1e-6);//
                double U = flux / inlet_lumen;
                /*        for (int I = 0; I < 10; I++)
                        {
                            inlet_lumen = Math.Pow(U - core_node.neighbours.Last().chrt, 4) / 64 * Math.Pow(Node.density / beta_1, 2);
                            U = flux / inlet_lumen;
                        }*/

                /*
                double inlet_lumen = myMath.NewtonSolver(flux_function, core_node.lumen_sq_0, 1e-10f, core_node.lumen_sq_0 * 1e-6);
                 * */

                double lambda = U + Math.Sqrt(beta_1 / 2 / Node.density) * Math.Pow(inlet_lumen, 0.25);
                double dx = dt * lambda;
                dx = dx * 10;

                double dZ = Vector3.Distance(core_node.neighbours.Last().position, core_node.position);
                double dUdt = (U - previous_velocity) / dt;
                double dUdZ = (core_node.neighbours.Last().velocity * v_sign[1] - U) / dZ;
                double visc_term = (double)(22 * Node.viscosity * Math.PI * U / Node.density / inlet_lumen);
                double P = core_node.neighbours.Last().pressure + (dUdt + U * dUdZ + visc_term) * dZ * Node.density;

                core_node.velocity = U * v_sign[0];
                core_node.pressure = P;
                core_node.lumen_sq = inlet_lumen;
                chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
                core_node.chrt = chrt;

                current_velocity = U;


                if (pulse_time_interval < Node.pulse_time_est)
                {
                    pulse_time_interval += dt;
                    pressure_hist.Enqueue(P);
                }
            }
        }       


        public               TableFunction flux_on_time;
        public readonly TableFunction base_flux_on_time;
        public readonly float               base_period;

        protected SimpleFunction flux_function;

        protected double beta_1;
        protected double Q_min;
        protected double Q_max;

        protected double previous_velocity;
        protected double current_velocity;
        protected double chrt;

        protected double diastolic_pressure;
        protected double sistolic_pressure;

        protected Queue<double> pressure_hist;

        protected double pulse_time_interval;
        protected double flux_minmax_time_interval;

        private float sample_dt;
    };


    public class PressureOutletRCR: BoundaryCondition
    {
        public PressureOutletRCR(BoundaryCondition BC, GetBetaFunction getElsticBeta, double _R1, double _R2, double _C)
            : base(BC.core_node, BC.current_time)
        {
            this.core_node = base.core_node;
            this.current_time = base.current_time;
            this.previous_time = base.previous_time;
            this.v_sign = base.v_sign;
            
            R1 = _R1;
            R2 = _R2;
             C =  _C;


            if (core_node.neighbours.Find(x => x.neighbours.Count <= 2) != null)
            {
                core_node.lumen_sq_0 = core_node.neighbours.Find(x => x.neighbours.Count <= 2).lumen_sq_0;
                core_node.lumen_sq   = core_node.lumen_sq_0;
                core_node.chrt       = core_node.neighbours.Find(x => x.neighbours.Count <= 2).chrt;
            }

            U_tx_10 = 0;
            U_tx_01 = 0;
            P_tx_10 = Node.diastolic_pressure;
            P_tx_01 = Node.diastolic_pressure;
            P_tx_11 = Node.diastolic_pressure;
            A_tx_01 = core_node.lumen_sq_0;
            A_tx_10 = core_node.lumen_sq_0;
            Q_t_1   = 0;
            Q_t_1   = 0;

            double R0 = Math.Sqrt(core_node.lumen_sq_0 / Math.PI);
            beta_1 = getElsticBeta(R0) / core_node.lumen_sq_0;
            flux_function = delegate(double A)
            {
                return (A * core_node.chrt - 4 * Math.Pow(A, 1.25f) * Math.Sqrt(beta_1 / 2.0 / Node.density) - Q_t_1) * 1000;
            };

            core_node.chrt = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
            chrt = core_node.chrt;
                      
            vel_function = delegate(double U_tx_11)
            {
                double dUdx = (U_tx_11 - U_tx_10) / interfunc_dx;
                double dUdt = (U_tx_11 - U_tx_01) / interfunc_dt;
                double dPdx = (P_tx_11 - P_tx_10) / interfunc_dx;
                double visc_term = (double)(22 * Node.viscosity * Math.PI * U_tx_10 /Node.density / core_node.neighbours.Last().lumen_sq);
                return dUdt + U_tx_11 * dUdx + 1 / Node.density * dPdx + visc_term;                
            };

            flux_function = delegate(double A_tx_11)
            {                
                double dAdt = (A_tx_11 - A_tx_01) / interfunc_dt;
                double dQdx = (A_tx_11 * U_tx_11 - A_tx_10 * U_tx_10) / interfunc_dx;
                return dAdt + dQdx;
            };

            DefineSign();
        }

        protected override void DefineSign()
        {
            Vector3 dir_vector1 = new Vector3();
            Vector3 dir_vector2 = new Vector3();


            dir_vector1 = core_node.position - core_node.neighbours.First().position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));

            dir_vector2 = core_node.neighbours.Last().dir_vector;
            v_sign[1] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));
        } 

        public override void doBC(double dt)
        {
            //if (false)
            {
                interfunc_dx = Vector3.Distance(core_node.position, core_node.neighbours.Last().position);
                interfunc_dt = dt;
                current_time = current_time + dt;
                previous_time = current_time;


                if (core_node.id == 1257)
                    core_node.id = 1257;

                Q_t_1 = core_node.lumen_sq * core_node.velocity * v_sign[0];

                double A_tx_11 = core_node.neighbours.Last().lumen_sq;//myMath.NewtonSolver(flux_function, core_node.lumen_sq_0, 1e-10f, core_node.lumen_sq_0 * 1e-6);

                double U_tx_11 = Q_t_1 / A_tx_11;

                //double dZ = Vector3.Distance(core_node.neighbours.Last().position, core_node.position);
                //double dUdt = (U_tx_11 - U_tx_01) / dt;
                //double dUdZ = (core_node.neighbours.Last().velocity * v_sign[1] - U_tx_11) / dZ;
                //double visc_term = (double)(22 * Node.viscosity * Math.PI * U_tx_11 / Node.density / A_tx_11);
                //double P_tx_11 = core_node.neighbours.Last().pressure - (dUdt + U_tx_11 * dUdZ + visc_term) * dZ * Node.density;
                        


                double     dQdt = (Q_t_1 - Q_t_0)  / dt;
                        P_tx_11 = (Q_t_0* (1 + R1 / R2) + C * R1 * dQdt + P_tx_01 * C / dt) / (C / dt + 1 / R2);

                

                core_node.velocity = U_tx_11 * v_sign[0];
                core_node.pressure = P_tx_11;
                core_node.lumen_sq = A_tx_11;
                chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);

                P_tx_01 = P_tx_11;
                A_tx_01 = A_tx_11;
                U_tx_01 = U_tx_11;
                Q_t_0 = Q_t_1;

                /*
                 double flux = flux_on_time(current_time);
                    double inlet_lumen = myMath.NewtonSolver(flux_function, core_node.lumen_sq_0, 1e-10f, core_node.lumen_sq_0*1e-6);
                    double U = flux / inlet_lumen;

                    double dZ = Vector3.Distance(core_node.neighbours.Last().position, core_node.position);           
                    double dUdt = (U - previous_velocity) / dt;
                    double dUdZ = (core_node.neighbours.Last().velocity*v_sign[1] - U) / dZ;
                    double visc_term = (double)(22 * Node.viscosity * Math.PI * U / Node.density / inlet_lumen);
                    double P = core_node.neighbours.Last().pressure + (dUdt + U * dUdZ + visc_term) * dZ * Node.density;
            
                    core_node.velocity = U*v_sign[0];
                    core_node.pressure = P;
                    core_node.lumen_sq = inlet_lumen;
                    chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);

                    current_velocity = U; 
                 //*/
            }
        }

        protected SimpleFunction flux_function;
        protected SimpleFunction  vel_function;
        
        protected double U_tx_01, U_tx_10, U_tx_11, P_tx_10, P_tx_11, P_tx_01;
        protected double A_tx_10, A_tx_01, A_tx_11, Q_t_0, Q_t_1;
        protected double interfunc_dx, interfunc_dt;

        public double R1, R2, C, beta_1;
        protected double chrt;

    }


    public class PressureOutletFree : BoundaryCondition
    {
        public PressureOutletFree(BoundaryCondition BC, GetBetaFunction getElsticBeta)
            : base(BC.core_node, BC.current_time)            
        {
            this.core_node = base.core_node;
            this.current_time = base.current_time;
            this.previous_time = base.previous_time;
            this.v_sign = base.v_sign;


            if (core_node.neighbours.Find(x => x.neighbours.Count <= 2) != null)
            {
                core_node.lumen_sq_0 = core_node.neighbours.Find(x => x.neighbours.Count <= 2).lumen_sq_0;
                core_node.lumen_sq = core_node.lumen_sq_0;
                core_node.chrt = core_node.neighbours.Find(x => x.neighbours.Count <= 2).chrt;
            }

            U_tx_10 = 0;
            U_tx_01 = 0;
            P_tx_10 = Node.diastolic_pressure;
            P_tx_01 = Node.diastolic_pressure;
            P_tx_11 = Node.diastolic_pressure;
            A_tx_01 = core_node.lumen_sq_0;
            A_tx_10 = core_node.lumen_sq_0;
            Q_t_1 = 0;
            Q_t_1 = 0;


            double R0 = Math.Sqrt(core_node.lumen_sq_0 / Math.PI);
            beta_1 = getElsticBeta(R0) / core_node.lumen_sq_0;
            flux_function = delegate(double A)
            {
                return (A * core_node.chrt - 4 * Math.Pow(A, 1.25f) * Math.Sqrt(beta_1 / 2.0 / Node.density) - Q_t_1) * 1000;
            };

            core_node.chrt = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);
            chrt = core_node.chrt;

            vel_function = delegate(double U_tx_11)
            {
                double dUdx = (U_tx_11 - U_tx_10) / interfunc_dx;
                double dUdt = (U_tx_11 - U_tx_01) / interfunc_dt;
                double dPdx = (P_tx_11 - P_tx_10) / interfunc_dx;
                double visc_term = (double)(22 * Node.viscosity * Math.PI * U_tx_10 / Node.density / core_node.neighbours.Last().lumen_sq);
                return dUdt + U_tx_11 * dUdx + 1 / Node.density * dPdx + visc_term;
            };

            flux_function = delegate(double A_tx_11)
            {
                double dAdt = (A_tx_11 - A_tx_01) / interfunc_dt;
                double dQdx = (A_tx_11 * U_tx_11 - A_tx_10 * U_tx_10) / interfunc_dx;
                return dAdt + dQdx;
            };

            DefineSign();
        }

        protected override void DefineSign()
        {
            Vector3 dir_vector1 = new Vector3();
            Vector3 dir_vector2 = new Vector3();


            dir_vector1 = core_node.position - core_node.neighbours.First().position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));

            dir_vector2 = core_node.neighbours.Last().dir_vector;
            v_sign[1] = Math.Sign(Vector3.Dot(dir_vector1, dir_vector2));
        }

        public override void doBC(double dt)
        {            
            {
                interfunc_dx = Vector3.Distance(core_node.position, core_node.neighbours.Last().position);
                interfunc_dt = dt;
                current_time = current_time + dt;
                previous_time = current_time;

                Q_t_1 = core_node.lumen_sq * core_node.velocity * v_sign[0];

                double A_tx_11 = core_node.neighbours.Last().lumen_sq;//myMath.NewtonSolver(flux_function, core_node.lumen_sq_0, 1e-10f, core_node.lumen_sq_0 * 1e-6);

                double U_tx_11 = Q_t_1 / A_tx_11;

                double dZ = Vector3.Distance(core_node.neighbours.Last().position, core_node.position);
                double dUdt = (U_tx_11 - U_tx_01) / dt;
                double dUdZ = (core_node.neighbours.Last().velocity * v_sign[1] - U_tx_11) / dZ;
                double visc_term = (double)(22 * Node.viscosity * Math.PI * U_tx_11 / Node.density / A_tx_11);
                double P_tx_11 = core_node.neighbours.Last().pressure;// -(dUdt + U_tx_11 * dUdZ + visc_term) * dZ * Node.density;

                core_node.velocity = U_tx_11 * v_sign[0];
                core_node.pressure = P_tx_11;
                core_node.lumen_sq = A_tx_11;
                chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / Node.density);

                P_tx_01 = P_tx_11;
                A_tx_01 = A_tx_11;
                U_tx_01 = U_tx_11;
                Q_t_0 = Q_t_1;
               
            }
        }

        protected SimpleFunction flux_function;
        protected SimpleFunction vel_function;

        protected double U_tx_01, U_tx_10, U_tx_11, P_tx_10, P_tx_11, P_tx_01;
        protected double A_tx_10, A_tx_01, A_tx_11, Q_t_0, Q_t_1;
        protected double interfunc_dx, interfunc_dt;

        protected double chrt, beta_1;
    }

    public class BC_spec
    {
        public BC_spec(BoundaryCondition bc_node)
        { 
            w_node = bc_node;
            Q_on_t = new List<double>();
            P_on_t = new List<double>();
              time = new List<double>();
     state_history = new List<PulseState>();
     CalcSpec();
        }

        public void Reset()
        {
            Q_on_t.Clear();
            P_on_t.Clear();
              time.Clear();
            period   = 0; 
        }

        public void Record(double _t)
        {
            if (w_node.core_node.id == 1257)
                w_node.core_node.id = 1257;
            Q_on_t.Add(-w_node.core_node.velocity * w_node.core_node.lumen_sq);
            P_on_t.Add(w_node.core_node.pressure);
            time.Add(_t);           
        }

        public void CalcSpec()
        {
            PulseState new_state = new PulseState ();
            try
            {
                new_state.Q_max = Q_on_t.Max();
                new_state.Q_min = Q_on_t.Min();

                new_state.Ps = P_on_t.Max();
                new_state.Pd = P_on_t.Min();
                int id = Q_on_t.FindIndex(x => x == new_state.Q_min);
                new_state.t_q_min = time[id];
                id = Q_on_t.FindIndex(x => x == new_state.Q_max);
                new_state.t_q_max = time[id];                

                new_state.td = time[P_on_t.FindIndex(x => x == new_state.Ps)];
                new_state.ts = time[P_on_t.FindIndex(x => x == new_state.Pd)];

                Q_tot = 0;
                period = 0;
                for (int i = 0; i < Q_on_t.Count() - 1; i++)
                {
                    Q_tot += (Q_on_t[i] + Q_on_t[i + 1]) / 2 * (time[i + 1] - time[i]);
                    period += (time[i + 1] - time[i]);
                }
                new_state.period = period;
                new_state.time_begin = time[0];
                new_state.time_end   = time[Q_on_t.Count() - 1];
                new_state.Q_average =  (Q_tot / period);
                new_state.Q_tot = Q_tot;

                this.state_history.Add(new_state);
            }
            catch (Exception e)
            {
                this.state_history.Add(new_state);
            }

            update();
        }

        public PulseState getState(int seq)
        {
            return state_history[state_history.Count - 1 - seq];
        }

        public void update()
        {
            Q_max = state_history.Last().Q_max;
            Q_min = state_history.Last().Q_min;
               Pd = state_history.Last().Pd;
               Ps = state_history.Last().Ps;

          t_q_max = state_history.Last().t_q_max;
          t_q_min = state_history.Last().t_q_min;

               td = state_history.Last().td;
               ts = state_history.Last().ts;
            Q_tot = state_history.Last().Q_tot;
         Q_averge = state_history.Last().Q_average;
           period = state_history.Last().period;
        }


        private BoundaryCondition w_node;

        public List<double> Q_on_t;
        public List<double> P_on_t;
        public List<double>   time;

        private List<PulseState> state_history;

        public double Q_max,     Q_min, Pd, Ps;
        public double t_q_max, t_q_min, td, ts;
        public double Q_tot,  Q_averge, period;       
    }

    public class RCRBCController
    {
        public RCRBCController(double _P_trg_dst, double _P_trg_sst, double _pulse_period)
        {
            P_trg_dst    = _P_trg_dst;
            P_trg_sst    = _P_trg_sst;
            P_mean = P_trg_dst + 1 / 3.0f * (P_trg_sst - P_trg_dst);

            pulse_period = _pulse_period;           

              start_time = 0;
               stop_time = 0;
             record_time = 0;
        }

        public void setPulsePeriod(double _pulse_period)
        {
            pulse_period = _pulse_period;
        }

        public void setBCset(List<PressureOutletRCR> _BC_set)
        {
            BC_set = _BC_set;
            BC_set_spec  = new List<BC_spec>();
            BC_spec_dict = new Dictionary<PressureOutletRCR, BC_spec>();
            foreach (var bc in _BC_set)
            {
                 BC_set_spec.Add(new BC_spec(bc));
                BC_spec_dict.Add(bc, BC_set_spec.Last());
            }
        }

        public void reset()
        {   
            foreach (var bc in BC_set_spec)
                bc.Reset();

            record_time = 0;
            start_time = 0;
            stop_time = 0;
        }
        
        public bool record(double curr_time)
        {
            if (start_time == 0)
                start_time = curr_time;

            if (record_time >= pulse_period - 1e-5)                            
                return false;            

            foreach (var bc in BC_set_spec)
                bc.Record(curr_time);

            record_time = curr_time - start_time;

            if (record_time >= pulse_period)
                stop_time = curr_time;            

            return true;
        }

        public bool getPressureConvergence(float relTol, out double av_SystPressure, out double av_DstPressure)
        {
            double max_pressure_diff  = 0;
            double max_pressure = 0;
            av_SystPressure = 0;
            av_DstPressure = 0;
            foreach (var bc in BC_set)
            {
                BC_spec bc_spec = BC_spec_dict[bc];
                bc_spec.CalcSpec();

                if(bc.core_node.id==1257)
                    bc_spec.CalcSpec();

                double prev_Pd = bc_spec.getState(1).Pd;
                double prev_Ps = bc_spec.getState(1).Ps;
                av_DstPressure += bc_spec.Pd;
                av_SystPressure += bc_spec.Ps;
                if (Math.Abs(prev_Pd - bc_spec.Pd) > max_pressure_diff)
                {
                    max_pressure_diff = Math.Abs(prev_Pd - bc_spec.Pd);
                    max_pressure = bc_spec.Ps;
                }
                if (Math.Abs(prev_Ps - bc_spec.Ps) > max_pressure_diff)
                {
                    max_pressure_diff = Math.Abs(prev_Ps - bc_spec.Ps);                    
                    max_pressure      = bc_spec.Ps;
                }
            }
            av_SystPressure = av_SystPressure / BC_set.Count;
            av_DstPressure = av_DstPressure / BC_set.Count;
            return (max_pressure_diff < max_pressure * relTol);                
        }

        public bool getFluxConvergence(float relTol)
        {
            double max_flux_diff = 0;
            double max_flux = 0;

            foreach (var bc in BC_set)
            {
                BC_spec bc_spec = BC_spec_dict[bc];
                bc_spec.CalcSpec();
                
                double prev_Q_max = bc_spec.getState(1).Q_max;
                double prev_Q_min = bc_spec.getState(1).Q_min;
                
                if (Math.Abs(prev_Q_min - bc_spec.Q_min) > max_flux_diff)
                {
                    max_flux_diff = Math.Abs(prev_Q_min - bc_spec.Q_min);
                    max_flux = Math.Abs(bc_spec.Q_max);
                }
                if (Math.Abs(prev_Q_max - bc_spec.Q_max) > max_flux_diff)
                {
                    max_flux_diff = Math.Abs(prev_Q_max - bc_spec.Q_max);
                    max_flux = Math.Abs(bc_spec.Q_max);
                }
            }
            return (max_flux_diff < relTol * max_flux);
        }

        public void adjustRCR()
        {
            foreach (var bc in BC_set)
            {
                BC_spec bc_spec = BC_spec_dict[bc];
                        bc_spec.CalcSpec();

                double c_d = Math.Sqrt(bc.beta_1 / 2.0f / Node.density) * Math.Pow(bc.core_node.lumen_sq_0, 0.25f);
                
                
                double P_mean = P_trg_dst + 1 / 3.0f * (P_trg_sst - P_trg_dst);

                double RT =  P_mean / bc_spec.Q_averge;                
                double R2 = 0;
                double R1 = 0;
                double C  = 0;

               // if (adj_iteration==1)
                {
                    C = (bc_spec.Q_max - bc_spec.Q_min) / (P_trg_sst - P_trg_dst) * Math.Abs(bc_spec.t_q_max - bc_spec.t_q_min);
                    R2 = RT - R1;
                    R1 = Node.density * c_d / bc.core_node.lumen_sq_0;
                    R2 = RT - R1;

                    if (R2 < 0)
                    {
                        continue;
                    }
                }
             /*   else
                {
                    C  = bc.C ; 
                    R1 = bc.R1;
                    R2 = bc.R2;
                }*/
                    
                
                int num_of_periods = 10;
                int dec = 3;
                int max_adj_cycles = 100;

                for (int N = 0; N < max_adj_cycles; N++)
                {
                    this.PressureRCRSimularot(num_of_periods, dec, bc_spec.Q_on_t, bc_spec.time, R1, R2, C);
                    double P_pulse = (simPs - simPd);
                    double P_pulse_trg = (P_trg_sst - P_trg_dst);

                    if (Math.Abs(P_pulse - P_pulse_trg) < 100 && Math.Abs(simPd - P_trg_dst) < 100)
                        break;

                    RT = RT + (P_trg_dst - simPd) / (bc_spec.Q_averge)*0.1;// *bc_spec.Q_averge / 10e-4;
                    R2 = RT - R1;
                    if (R2 < 0)
                    {
                        break;                        
                    }

                    C = C + (bc_spec.Q_max - bc_spec.Q_min) / (P_pulse * P_pulse) * Math.Abs(bc_spec.t_q_max - bc_spec.t_q_min) * (P_pulse - P_pulse_trg) * 0.1;
                    if (C < 0)                  
                        C = 1e-15;                 
                  
                }              
                
                    bc.C = C;
                    bc.R1 = R1;
                    bc.R2 = R2;                
            }
        }

        public void PressureRCRSimularot(int num_of_periods, int decimation, List<double> Q_on_t, List<double> time, 
                                         double R1, double R2, double C)
        {
            double P_next = P_trg_dst; double P_curr = P_trg_dst;

            simPd = double.MaxValue;
            simPs = 0;

            for (int N = 0; N < num_of_periods; N++)
            {
                for (int i = decimation; i < Q_on_t.Count - decimation; i += decimation)
                {
                    P_curr = P_next;
                    double dQdt = (Q_on_t[i + decimation] - Q_on_t[i - decimation]) / (time[i + decimation] - time[i - decimation]);
                    double dt = (time[i + decimation] - time[i]);
                    P_next = ((Q_on_t[i] * (R2 + R1) + C * R1 * R2 * dQdt) * dt + C * R2 * P_curr) / (dt + C * R2);

                  //  P_next = (Q_on_t[i] * (1 + R1 / R2) + C * R1 * dQdt + P_curr * C / dt) / (C / dt + 1 / R2);

                    if (N == num_of_periods - 1)
                    {
                        if (P_curr > simPs)
                        {
                            simPs = P_curr;
                            sim_ts = time[i];
                        }

                        if (P_curr < simPd)
                        {
                            simPd = P_curr;
                            sim_td = time[i];
                        }
                    }
                }
            }
        }



        public double  record_time;
        public double   start_time;
        public double    stop_time;


        private List <PressureOutletRCR> BC_set;
        private List <BC_spec>           BC_set_spec;        
        private Dictionary<PressureOutletRCR,BC_spec> BC_spec_dict;
        private double P_trg_dst, P_trg_sst, P_mean, pulse_period;
        private double simPd, simPs, sim_td, sim_ts;      
    }

    public delegate int intSystemMask(Node n);

    public struct BC_params
    {
        public int id;
        public double R1, R2, C;
    }

    public class IO_Module
    {
        static public bool LoadTableFunction(string filename, double period, out TableFunction t_f)
		{
            return LoadTableFunctionFromString( File.ReadAllText(filename), period, out t_f );
		}


        static public bool LoadTableFunctionFromString(string text, double period, out TableFunction t_f)
        {
            t_f = null;

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo(localization);

            string[] readText = text.Split(new[]{'\r','\n'}, StringSplitOptions.RemoveEmptyEntries);
            if (readText.GetLength(0) == 0)
                return false;

            bool res = false;
            double[,] table_function_0 = new double[readText.Length, 2];
            double[,] table_function_1 = new double[readText.Length, 2];
            Regex regex = new Regex(@"^\s*(-*\d+(\.\d+)?(e-*\d+)?)$", RegexOptions.IgnoreCase);
            Match node_match = regex.Match(readText[0]);
            double scale = double.Parse(node_match.Groups[1].Value);
            for (int i = 1; i < table_function_1.GetLength(0); i++)
            {
                regex = new Regex(@"^\s*(-*\d+\.\d+)[\s+\t]?(-*\d+\.\d+)$", RegexOptions.IgnoreCase);
                node_match = regex.Match(readText[i]);
                double time = double.Parse(node_match.Groups[1].Value);
                double val  = double.Parse(node_match.Groups[2].Value);
                table_function_1[i, 0] = time;
                table_function_1[i, 1] = val;
                table_function_0[i, 0] = time;
                table_function_0[i, 1] = val;
                res = true;
            };

            if (res)
            {
                t_f = delegate(double time)
                {
                    double value = 0;

                    int l_ind = 0;
                    int r_ind = table_function_0.GetLength(0) - 1;

                    time = time - period * Math.Floor(time / period);
                    {

                        while ((r_ind - l_ind) > 1)
                        {
                            int cur_ind = l_ind + (r_ind - l_ind) / 2;
                            if (table_function_0[cur_ind, 0] > time)
                                r_ind = cur_ind;
                            else
                            {
                                l_ind = cur_ind;
                            }
                        }

                        double l_time = table_function_0[l_ind, 0];
                        double r_time = table_function_0[r_ind, 0];

                        double l_value = table_function_0[l_ind, 1];
                        double r_value = table_function_0[r_ind, 1];

                        value = l_value + (time - l_time) / (r_time - l_time) * (r_value - l_value);
                    }

                    return value * scale;
                };
            }

            return res;
        }

        static public TableFunction makeTableFunction(double[,] table_function_0)
        {
                double period = table_function_0[table_function_0.GetLength(0) - 1, 0] - table_function_0[0, 0];

                TableFunction t_f = delegate(double time)
                {
                    double value = 0;

                    int l_ind = 0;
                    int r_ind = table_function_0.GetLength(0) - 1;

                    time = time - period * Math.Floor(time / period);
                    {

                        while ((r_ind - l_ind) > 1)
                        {
                            int cur_ind = l_ind + (r_ind - l_ind) / 2;
                            if (table_function_0[cur_ind, 0] > time)
                                r_ind = cur_ind;
                            else
                            {
                                l_ind = cur_ind;
                            }
                        }

                        double l_time = table_function_0[l_ind, 0];
                        double r_time = table_function_0[r_ind, 0];

                        double l_value = table_function_0[l_ind, 1];
                        double r_value = table_function_0[r_ind, 1];

                        value = l_value + (time - l_time) / (r_time - l_time) * (r_value - l_value);
                    }

                    return value;
             };
             return t_f;
        }


        static public TableFunction xScaleTableFunction(double timestep, double period, double new_period, TableFunction t_f)
        {   
            int N = (int)Math.Ceiling(period/timestep);
            timestep = period / N;
            double[,] table_function = new double[N, 2];

            double scale = new_period / period;

            for (int i = 0; i < N; i++)
            {
                table_function[i, 0] = timestep * i * scale;
                table_function[i, 1] = t_f(timestep * i) / scale;
            }
            return IO_Module.makeTableFunction(table_function);
        }

        static public void WriteRCR(string filename, List <PressureOutletRCR> bc_list)
        {
            List<string> writeText = new List<string>();
            for (int i = 0; i < bc_list.Count; i++)
            {
                string out_text = "";
                out_text = out_text + bc_list[i].core_node.id.ToString() + " R1:" + (bc_list[i].R1*1e-9).ToString() + " R2:" + (bc_list[i].R2*1e-9).ToString() + " C:" + (bc_list[i].C*1e12).ToString();
                writeText.Add(out_text);
            }
            File.WriteAllLines(filename, writeText);
        }

        static public bool LoadBC_RCR_params(string filename, out List <BC_params> BC_Params)
		{
			return LoadBC_RCR_paramsFromString( File.ReadAllText(filename), out BC_Params );
		}

        static public bool LoadBC_RCR_paramsFromString(string text, out List <BC_params> BC_Params)
        {
            BC_Params = new List<BC_params>();

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo(localization);
            string[] readText = text.Split(new[]{'\r','\n'}, StringSplitOptions.RemoveEmptyEntries);

            bool res = false;
            Match node_match;
            for (int i = 0; i < readText.GetLength(0); i++)
            {
                                            //56          R1:  1.8104          R2:  7.2417         C:   31.29
                Regex regex = new Regex(@"^(\d+)\s+R1:(\d+\.\d+)\s+R2:(\d+\.\d+)\s+C:(\d+\.?\d*)$", RegexOptions.IgnoreCase);//^(\d+)[\s+\t]?R1:(\d+\.\d+)[\s+\t]?R2:(\d+\.\d+)[\s+\t]?C:(\d+\.\d+)$", RegexOptions.IgnoreCase);
                node_match = regex.Match(readText[i]);
                int      id =    int.Parse(node_match.Groups[1].Value);
                double   R1 = double.Parse(node_match.Groups[2].Value);
                double   R2 = double.Parse(node_match.Groups[3].Value);
                double    C = double.Parse(node_match.Groups[4].Value);

                BC_params prms = new BC_params();
                BC_Params.Add(prms);

                prms.id = id;
                prms.R1 = R1*1e9;
                prms.R2 = R2*1e9;
                prms.C  = C*1E-12;

                BC_Params[i] = prms;

                res = true;
            };

            return res;
        }

        public static string localization;
    }

    public class VascularNet
    {
        static public int NumOfNeigbours(Node node)
        {
            return node.neighbours.Count;
        }

        public double setHeartRate(int HR, double timestep)
        {
            if (HR > 120 || HR < 40)
                return 0;

            double new_period = 60.0/HR;

            foreach(var bc in this.bounds)
                if (bc.GetType() == typeof(InletFlux))
                {
                    InletFlux inlet_bc = (InletFlux)bc;
                    inlet_bc.flux_on_time = IO_Module.xScaleTableFunction(timestep, inlet_bc.base_period, new_period, inlet_bc.base_flux_on_time);
                }

            return new_period;
        }


        public static string LoadFromFile(VascularNet vnet, String file_path)
		{
			return LoadFromString( vnet, File.ReadAllText(file_path) );
		}			 


		public static string LoadFromString (VascularNet vnet, String text)
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-US");
            vnet.vascular_system.Clear();

            string protocol = "VescularNet loading protocol:\n";

            string[] readText = text.Split(new[]{'\r','\n'});

            Regex regex = new Regex(@"^name:\s*(\w+)$", RegexOptions.IgnoreCase);
            int i = 0;
            while (!regex.IsMatch(readText[i]))
            {
                i++;
                if (i >= readText.Length)
                {
                    protocol += "Error: No correct name string was found!\n";
                    return protocol;
                }
            }

            Match name_match = regex.Match(readText[i]);
            vnet.name = name_match.Groups[1].Value;

            protocol += "The name was read: " + vnet.name;
            protocol += ";\n";

            Regex regex_1 = new Regex(@"^Coordinates:\s*$", RegexOptions.IgnoreCase);
            Regex regex_2 = new Regex(@"^Bonds:\s*$", RegexOptions.IgnoreCase);

            List<List<int>> bonds_index = new List<List<int>>();
            int node_count = 0;
            int bond_string_count = 0;

            while (true)
            {

                i++;

                if (regex_1.IsMatch(readText[i]))
                    while (true)
                    {
                        i++;
                        if (i >= readText.Length)
                            break;


                        regex = new Regex(@"^\s*(\d+)\s+X:(-*\d+.\d+)\s+Y:(-*\d+.\d+)\s+Z:(-*\d+.\d+)\s+R:(-*\d+.\d+)\s+C:(\d+.\d+)$", RegexOptions.IgnoreCase);
                        regex = new Regex(@"^\s*(\d+)\s+X:(-*\d+.\d+)\s+Y:(-*\d+.\d+)\s+Z:(-*\d+.\d+)\s+R:(-*\d+.\d+)\s+C:(\d+.\d+)$", RegexOptions.IgnoreCase);
                        Match node_match = regex.Match(readText[i]);

                        if (node_match.Groups.Count < 6)
                        {
                            regex = new Regex(@"^\s*(\d+)\s+X:(-*\d+.\d+)\s+Y:(-*\d+.\d+)\s+Z:(-*\d+.\d+)\s+R:(-*\d+.\d+)$", RegexOptions.IgnoreCase);
                            node_match = regex.Match(readText[i]);
                            if (node_match.Groups.Count < 6)
                                break;
                        }


                        int id = int.Parse(node_match.Groups[1].Value);
                        Vector3 position = new Vector3(double.Parse(node_match.Groups[2].Value),
                                                       double.Parse(node_match.Groups[3].Value),
                                                       double.Parse(node_match.Groups[4].Value));
                        double rad = double.Parse(node_match.Groups[5].Value);

                        vnet.vascular_system.Add(new Node(id, position, rad * 1e-3));
                        node_count++;
                    }
                else
                    if (regex_2.IsMatch(readText[i]))
                        while (true)
                        {
                            i++;
                            if (i >= readText.Length)
                                break;

                            regex = new Regex(@"(\d+)\s*", RegexOptions.IgnoreCase);
                            MatchCollection node_match = regex.Matches(readText[i]);
                            if (node_match.Count < 2)
                                break;


                            int id = int.Parse(node_match[0].Value);
                            bonds_index.Add(new List<int>());
                            bonds_index[bonds_index.Count - 1].Add(id);

                            for (int n = 1; n < node_match.Count; n++)
                                bonds_index[bonds_index.Count - 1].Add(int.Parse(node_match[n].Value));

                            bond_string_count++;
                        }

                if (i >= readText.Length - 1)
                {
                    protocol += node_count.ToString() + " nodes were read;\n";
                    protocol += bond_string_count.ToString() + " bonds strings were read;\n";

                    foreach (var str in bonds_index)
                    {
                        Node nd = vnet.vascular_system.Find(x => x.id == str[0]);
                        for (int s = 1; s < str.Count; s++)
                            nd.addNeighbours(new Node[] { vnet.vascular_system.Find(x => x.id == str[s]) });
                    }

                    break;
                }
            }

            return protocol;
        }

        public VascularNet(string _name)
        {
            name = new string(_name.ToCharArray());
            vascular_system = new List<Node>();

            knots = new List<Knot>();
            bounds = new List<BoundaryCondition>();
            threads = new List<Thread>();
        }

        public void setCloth(int node_id, float degree)
        {
            Node    nd = this.vascular_system.Find(x => x.id == node_id);
            Thread thr = this.node2thread[nd];

            thr.setClot(nd, degree);           
        }

        public void removeCloth(int node_id)
        {
            Node nd = this.vascular_system.Find(x => x.id == node_id);
            Thread thr = this.node2thread[nd];
            thr.removeClot(nd);           
        }

        public void WriteThread(Thread thread)
        {
            string out_text = "";
            out_text += "Name: System_0\n";
            out_text += "Coordinates:\n";
            foreach (var n in thread.nodes)
            {
                out_text += n.id.ToString();
                out_text += " X:" + n.position.x.ToString("F5") + " Y:" + n.position.y.ToString("F5") + " Z:" + n.position.z.ToString("F5");
                out_text += " R:" + (Math.Sqrt(n.lumen_sq_0 / Math.PI) * 1e+3).ToString("F5") + "\n";
            }
            out_text += "Bonds:\n";
            foreach (var n in thread.nodes)
            {
                out_text += n.id.ToString();
                foreach (var nn in n.neighbours)
                    out_text += " " + nn.id.ToString();
                out_text += "\n";
            }
            System.IO.File.WriteAllText(@"thread.txt", out_text);
        }

        public List<Node> getSubsystem(NodeFilter filter, int val)
        {

            if (filter == null)
                return vascular_system.GetRange(0, vascular_system.Count);
            else
            {
                List<Node> subsystem = new List<Node>();

                foreach (Node n in vascular_system)
                    if (filter(n) == val)
                        subsystem.Add(n);
                return subsystem;
            }
            return null;
        }

        public List<Node> getSubsystem(Predicate<Node> p)
        {

            if (p == null)
                return vascular_system.GetRange(0, vascular_system.Count);
            else
            {
                List<Node> subsystem = new List<Node>();

                foreach (Node n in vascular_system)
                    if (p(n))
                        subsystem.Add(n);
                return subsystem;
            }
            return null;
        }

        public void setSubsystem(List<Node> sub_sustem)
        {
            foreach (var n in sub_sustem)
            {
                int ind = vascular_system.FindIndex(x => x.id == n.id);
                List<Node> neughbours = vascular_system[ind].neighbours;
                n.addNeighbours(neughbours.ToArray());
                foreach (var nn in neughbours)
                {
                    nn.addNeighbour(n);
                    nn.neighbours.Remove(vascular_system[ind]);
                }

                vascular_system[ind] = n;
            }
        }
        
        public void defineNet()
        {
            int count_id = vascular_system.Count;

            foreach (var n in vascular_system)
                n.defDirVector();

            List<Node> knot_nodes = getSubsystem(x => x.neighbours.Count >  2);
            List<Node> term_nodes = getSubsystem(x => x.neighbours.Count == 1);

            knots  = new List<Knot>();
            bounds = new List<BoundaryCondition>();

            foreach (var kn in knot_nodes)
            {
                List<Node> b_nodes = new List<Node>();
                List<Node> n_nodes = new List<Node>(kn.neighbours);

                foreach (var n in n_nodes)
                {
                    Vector3 pos = kn.position - (kn.position - n.position) * 0.01;
                    double    r = Math.Sqrt(n.lumen_sq_0/Math.PI);
                    Node b_n = new Node(count_id, pos, r);
                    b_n.neighbours.Add(kn); b_n.neighbours.Add(n);
                    kn.neighbours.Remove(n);
                    n.neighbours.Remove(kn);
                    kn.addNeighbour(b_n);
                    n.addNeighbour(b_n);
                    b_nodes.Add(b_n);
                    count_id++;
                }
                knots.Add(new Knot(kn, 0));
            }

            foreach(var tn in term_nodes)
                bounds.Add(new BoundaryCondition(tn, 0));

            threads     = new List<Thread>();
            node2thread = new Dictionary<Node, Thread>();

            getBoolValueDelegate isProcessed;
            setBoolValueDelegate setProcessed;
            Node.newBoolValueLayer(out isProcessed, out setProcessed);

            Node curr_node;
            Node next_node;

            foreach (var kn in knots)
                setProcessed(kn.core_node, true);

            foreach (var bc in bounds)
                setProcessed(bc.core_node, true);
            
            foreach(var kn in knots)
            {   
                foreach (var n in kn.nodes)
                {
                    List<Node> protothread = new List<Node>();
                    curr_node = kn.core_node;
                 //   protothread.Add(curr_node);

                    next_node = n;
                    if (isProcessed(next_node) && next_node.neighbours.Count > 1)
                        continue;                    
                    
                    while(true)
                    {   
                        setProcessed   (curr_node, true);
                        protothread.Add(curr_node);
                        next_node = curr_node.neighbours.Find(x => !isProcessed(x));

                        if (next_node !=null&& next_node.id == 23382)
                            setProcessed(curr_node, true);

                        if (next_node == null)
                        {
                            next_node = curr_node.neighbours.Find(x => x.neighbours.Count != 2);
                            if (next_node != null)
                            {
                                protothread.Add(next_node);
                                if (protothread[0].neighbours.Count>2)
                                    protothread.RemoveAt(0);
                                if (protothread[protothread.Count-1].neighbours.Count > 2)
                                    protothread.RemoveAt(protothread.Count - 1);

                                Thread thread = new Thread(protothread);                               
                                threads.Add(thread);                              
                            }
                            break;
                        }
                        curr_node = next_node;
                    }
                }
            }

            foreach (var tr in threads)
            {
                int WIN_LENGTH = 5;
                double[] window = new double[WIN_LENGTH];

                if (tr.nodes.GetLength(0) < WIN_LENGTH)
                    continue;

                for (int i = 0; i < tr.nodes.GetLength(0); i++)
                {
                    int offset = WIN_LENGTH / 2;
                    if (i < WIN_LENGTH / 2)
                        offset = i;

                    if (i >= tr.nodes.GetLength(0) - WIN_LENGTH / 2 )
                        offset = WIN_LENGTH - (tr.nodes.GetLength(0) - i);
                    
                    for (int j = 0; j < WIN_LENGTH; j++)
                        window[j] = tr.nodes[i - offset + j].lumen_sq_0;

                    window = window.OrderBy(x => x).ToArray();
                    tr.nodes[i].lumen_sq_0 = window[WIN_LENGTH / 2];
                    tr.nodes[i].lumen_sq   = window[WIN_LENGTH / 2];
                }
            }

            foreach (var bc in bounds)
            {
                List<Node> protothread = new List<Node>();
                
                curr_node = bc.core_node;
                protothread.Add(curr_node);

                next_node = curr_node.neighbours.Last();
                if (isProcessed(next_node))
                    continue;

                while (true)
                {  
                    curr_node = next_node;
                    setProcessed(curr_node, true);
                    protothread.Add(curr_node);
                    next_node = curr_node.neighbours.Find(x => !isProcessed(x));

                    if (next_node == null)
                    {
                        next_node = curr_node.neighbours.Find(x => x.neighbours.Count != 2);
                        if (next_node != null)
                        {
                            protothread.Add(next_node);
                            Thread thread = new Thread(protothread);                            
                            threads.Add(thread);
                        }
                        break;
                    }
                }
            }



            Node.terminateBoolValueLayer(ref isProcessed);
        }

        public List<BoundaryCondition>  getBounds()
        {
            return bounds;
        }

        public List<Thread> getThreads()
        {
            return threads;
        }

        public void holdThread()
        {                   
            foreach (var tr in threads)
                foreach (var nd in tr.nodes)
                    this.node2thread.Add(nd, tr);
        }

        private string name;
        private List<Node>      vascular_system;

        public List<Thread>             threads;
        public List<Knot>                 knots;
        public List<BoundaryCondition>   bounds;

        public Dictionary<Node, Thread> node2thread;


    };

}


