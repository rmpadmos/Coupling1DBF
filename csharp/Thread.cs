#define SAFE
//#define FAST

using System;
using System.Collections.Generic;
using System.Linq;

namespace BloodFlow
{

    public struct Clot
    {
        public VascularNode c_nd;
        public List<VascularNode> nd;
        public int thread_id;
        public double tgr_degree;
        public double curr_degree;
        public List<double> normal_lumen;
    }

    public class FFRClot       
    {
        public FFRClot(double curr_time, double _alpha)
        {
            p_pressure_acc = 0;
            d_pressure_acc = 0;
            av_flux_acc = 0;
            last_time = curr_time;
            alpha = _alpha;
        }

        public void init()
        {
            p_pressure_acc = 0;
            d_pressure_acc = 0;
            av_flux_acc = 0;
            
            length = 0;
            for (int i = 1; i < nd.Count; i++)
                length += Vector1.Length(nd[i].position - nd[i - 1].position);

            effective_visc = GlobalDefs.BLOOD_VISC;
            tg_ffr = 0.5;
        }

        public bool update_clot(double curr_time, double heart_period)
        {
            double dt = curr_time - last_time;
            last_time = curr_time;

            p_pressure_acc += proximal_nd.pressure * dt;
            d_pressure_acc +=   distal_nd.pressure * dt;
            av_flux_acc += proximal_nd.velocity * proximal_nd.lumen_sq * dt;

            if (curr_time % heart_period <= dt)
            {
                proximal_pressure = p_pressure_acc / heart_period;
                p_pressure_acc = 0;
                  distal_pressure = d_pressure_acc / heart_period;
                d_pressure_acc = 0;

                av_flux = Math.Abs(av_flux_acc);                
                av_flux_acc = 0;

                  if (proximal_pressure < distal_pressure)
                  {
                      VascularNode tmp = proximal_nd;
                      double tmp_p = proximal_pressure;
                      proximal_nd = distal_nd;
                      distal_nd = tmp;
                      proximal_pressure = distal_pressure;
                      distal_pressure = tmp_p;                      
                  }

                curr_ffr = distal_pressure / proximal_pressure;
                effective_visc = effective_visc * (1 + 4 * (curr_ffr - tg_ffr));
                
                //(1 - tg_ffr) / length / alpha / Math.PI * proximal_pressure / av_flux;
                

                return true;
            }

            return false;
        }

        public double   tg_ffr;
        public double curr_ffr;
        public double effective_visc;

        public List<VascularNode> nd;
        public VascularNode proximal_nd;
        public VascularNode   distal_nd;
        public int thread_id;
        public double proximal_pressure;
        public double   distal_pressure;        
        public double av_flux;

        protected double p_pressure_acc;
        protected double d_pressure_acc;
        protected double av_flux_acc;
        protected double length;
        protected double alpha;

        protected double last_time;
    }

    public class Thread
    {
        public Thread(List<VascularNode> protothread)
        {


            int L = protothread.Count;
            nodes = (VascularNode[])protothread.ToArray().Clone();
            v_sign = new int[L];

            velocity     = new double[L];
            pressure     = new double[L];
            lumen_sq     = new double[L];
            flux = new double[L];
            //chrt = new double[L];

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
            }
        }

        public virtual void updateState()
        { }

        public virtual void updateStateFFR()
        { }

        unsafe public void state2nodes()
        {
            int L = nodes.GetLength(0);
            for (int i = 0; i < L; i++)
            {
                nodes[i].velocity = velocity[i] * v_sign[i];
                nodes[i].pressure = pressure[i];
                nodes[i].lumen_sq = lumen_sq[i];             
            }
        }

        public int getLength()
        {
            return nodes.GetLength(0);
        }

        public virtual void calcThread(double dt)
        {
        }

        public virtual void reset()
        {
           
        }

        protected void DefineSigns()
        {
            int L = nodes.GetLength(0);  

            Vector1[] dir_vector1 = new Vector1[L];
            Vector1[] dir_vector2 = new Vector1[L];
            v_sign = new int[L];
            for (int i = 0; i < L - 1; i++)
            {
                dir_vector1[i] = (nodes[i + 1].position - nodes[i].position);
            }
            for (int i = 0; i < L; i++)
            {
                nodes[i].defDirVector();
                dir_vector2[i] = nodes[i].dir_vector;
            }

            if (L > 1)
                dir_vector1[L - 1] = nodes[L - 1].position - nodes[L - 2].position;
            else
                dir_vector1[L - 1] = dir_vector2[L - 1];

            for (int i = 0; i < L; i++)
                v_sign[i] = Math.Sign(Vector1.Dot(dir_vector1[i], dir_vector2[i]));
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
                outText += n.id + " X:" + n.position.x.ToString("F8") + " R:" + Math.Sqrt(n.lumen_sq_0 / Math.PI).ToString("F8") + "\n";
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

        public virtual bool setClot(VascularNode nd, float degree)
        { return false; }

        public virtual bool setFFRClot(VascularNode nd, float FFR_degree)
        { return false; }

        public virtual void removeClot(VascularNode nd)
        { }

     /*   public void clearClotList()
        {
            clotes.Clear();
        }*/


        public VascularNode[] nodes;
        public int[] v_sign;

        public double[] velocity;
        public double[] pressure;
        public double[] lumen_sq;
        public double[] flux;
        
        public double[] chrt;

        public double current_time;
        protected List<Clot>     clotes;
        protected List<FFRClot>  ffr_clotes;
    }


    public class ElasticThread : Thread
    {
        public ElasticThread(Thread _protothread, GetBetaFunction getElsticBeta)
            : base(_protothread.nodes.ToList())
        {
         

            int L = nodes.GetLength(0);
            beta_1 = new double[L];
            lumen_sq_0 = new double[L];
            wall_thickhess = new double[L];
            viscosity = new double[L];            
            g_energy = new double[L];


            for (int i = 0; i < L; i++)
            {

                double R0 = Math.Sqrt(nodes[i].lumen_sq_0 / Math.PI);
                beta = getElsticBeta(R0, nodes[i].elasticity);

                beta_1[i] = beta / nodes[i].lumen_sq_0;
                wall_thickhess[i] = GlobalDefs.getBoileauWallThickness(R0, nodes[i].elasticity);
                lumen_sq_0[i] = nodes[i].lumen_sq_0;
                lumen_sq[i] = nodes[i].lumen_sq_0;
                pressure[i] = GlobalDefs.DIASTOLIC_PRESSURE;                
                //g_energy[i] = Vector1.Dot(nodes[i].position - GlobalDefs.ZERO_POINT, GlobalDefs.DOWN) * GlobalDefs.GRAVITY;
                g_energy[i] = 0;
                viscosity[i] = GlobalDefs.BLOOD_VISC;
            }
#if FAST
            unsafe
            {
                double* dZ = (double*)Marshal.AllocHGlobal((L - 1) * sizeof(double));
                for (int i = 0; i < L - 1; i++)                
                    dZ[i] = Vector3.Distance(nodes[i + 1].position, nodes[i].position);


                thread_solver = new McCormackThread(L, dZ, GlobalDefs.BLOOD_DENSITY, GlobalDefs.BLOOD_VISC, GlobalDefs.FRICTION_C);

                for (int i = 0; i < L; i++)
                {                    
                    int I = i;

                    _1DFunction_del pressure_func_del = delegate(double x)
                    {
                        return calcPressure(x, I);
                    };

                    delegate1DFunc f = new delegate1DFunc(pressure_func_del);

                    thread_solver.addFunc(f, i);

                }

                fixed (double* g_ptr = &g_energy[0])
                {
                    //thread_solver.setGravity(g_ptr);
                }

            }

#endif

            clotes = new List<Clot>();
            ffr_clotes = new List<FFRClot>();   
        }

        public override void reset()
        {
            int L = nodes.GetLength(0);

            foreach (var clt in clotes)
            {
                int CL = 2; int sh = 0;
                int nd_id = clt.thread_id;
                for (int j = -CL; j < CL; j++)
                    if (nd_id + j >= 0 && nd_id + j < nodes.GetLength(0))
                    {
                        clt.nd[j + CL + sh].lumen_sq_0 = clt.normal_lumen[j + CL + sh];
                        lumen_sq_0[nd_id + j] = clt.nd[j + CL + sh].lumen_sq_0;
                    }
                    else
                        sh--;
            }
                

            for (int i = 0; i < L; i++)
            {
                velocity[i] = 0;
                pressure[i] = 0;                
                lumen_sq[i] = lumen_sq_0[i];
                
            }

            clotes.Clear();

            state2nodes();
        }


        public override bool setClot(VascularNode nd, float degree)
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
                clt.nd = new List<VascularNode>();

                int L = 2;

                for (int i = -L; i < L+1; i++)
                    if (nd_id + i >= 0 && nd_id + i < nodes.GetLength(0))
                    {
                        clt.nd.Add(nodes[nd_id + i]);
                        clt.normal_lumen.Add(lumen_sq_0[nd_id + i]);
                    }

                clt.tgr_degree = degree;
                clt.curr_degree = 0;
                clt.thread_id = nd_id;
                clotes.Add(clt);
                return true;
            }

            return false;
        }


        public override bool setFFRClot(VascularNode nd, float FFR_degree)
        {
            int nd_id = 0;
            for (int i = 0; i < nodes.GetLength(0); i++)
                if (nodes[i] == nd)
                { nd_id = i; break; }

            if (clotes.FindIndex(x => x.c_nd == nd) == -1)
            {

                FFRClot clt = new FFRClot(current_time, GlobalDefs.FRICTION_C);
                
                clt.nd = new List<VascularNode>();
                int L = 2;
                if (nd_id - L - 1 >= 0)
                    clt.proximal_nd = nodes[nd_id - L - 1];
                else
                    return false;

                if (nd_id + L + 1 < nodes.GetLength(0))
                    clt.distal_nd = nodes[nd_id + L + 1];
                else
                    return false;                

                for (int i = -L; i < L + 1; i++)
                    if (nd_id + i >= 0 && nd_id + i < nodes.GetLength(0))
                    {                       
                        clt.nd.Add(nodes[nd_id + i]);                       
                    }

                clt.thread_id = nd_id;
                clt.init();
                ffr_clotes.Add(clt);
                return true;
            }

            return false;
        }


        public override void removeClot(VascularNode nd)
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
                if (Math.Abs(clt.tgr_degree - clt.curr_degree) != 0.0)
                {
                    clt.curr_degree += Math.Sign(clt.tgr_degree - clt.curr_degree)*0.001f;

                    if (Math.Abs(clt.tgr_degree - clt.curr_degree) < 0.01)
                        clt.curr_degree = clt.tgr_degree;

                    int L = 2; int sh = 0;
                    int nd_id = clt.thread_id;
                    for (int j = -L; j < L+1; j++)
                    {
                        if (nd_id + j >= 0 && nd_id + j < nodes.GetLength(0))
                        {
                            double degree = clt.curr_degree * (L - Math.Abs(j) + 1) / (L + 1);
                            clt.nd[j + L + sh].lumen_sq_0 = clt.normal_lumen[j + L + sh] * (1 - degree);                            
                            lumen_sq_0[nd_id + j] = clt.nd[j + L + sh].lumen_sq_0;
                            beta_1[nd_id + j] = this.beta / clt.nd[j + L + sh].lumen_sq_0;
                        }
                        else
                            sh--;
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
            curr_dt = dt;
            int L = nodes.GetLength(0);

            this.nodes2state();

            //if (nodes.Last().id == 13234 || nodes.First().id == 13234)
            //{
            //    L = nodes.GetLength(0);
            //}

        //    lumen_sq_old = (double[])lumen_sq.Clone();
#if FAST
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
#endif

      //     this.nodes2state();
#if SAFE       
          
            
           double[] velocity_pred = (double[])velocity.Clone();
           double[] lumen_sq_pred = (double[])lumen_sq.Clone();
           double[] pressure_pred = (double[])pressure.Clone();

           for (int i = 1; i < L-1; i++)
           {
               dz = Vector1.Distance(nodes[i].position, nodes[i - 1].position);
				//dz = 2.5e-3;
				velocity_pred[i] = velocity[i] - dt / dz * ((velocity[i] * velocity[i] / 2 + pressure[i] / GlobalDefs.BLOOD_DENSITY + g_energy[i]) - (velocity[i - 1] * velocity[i - 1] / 2 + pressure[i - 1] / GlobalDefs.BLOOD_DENSITY + g_energy[i-1]));
               velocity_pred[i] = velocity_pred[i] - 1.0 / GlobalDefs.BLOOD_DENSITY / lumen_sq[i] * (GlobalDefs.FRICTION_C * viscosity[i] * Math.PI * velocity_pred[i]) * dt;
               lumen_sq_pred[i] = lumen_sq[i] - dt / dz * (lumen_sq[i] * velocity[i] - lumen_sq[i-1] * velocity[i-1]);
               pressure_pred[i] = calcPressure(lumen_sq_pred[i],i);
                //Console.WriteLine(i);
                   
           }           
           //for(int i=0; i<L;i++)
            //{
            //    if (nodes[i].blocked == true)
            //        velocity_pred[i] = 0;
            //}

            for (int i = 1; i < L-1; i++)
           {
			    dz = Vector1.Distance(nodes[i].position, nodes[i + 1].position);
				//dz = 2.5e-3;
               velocity[i] = (velocity[i] + velocity_pred[i]) / 2 - dt / dz / 2 * ((velocity_pred[i + 1] * velocity_pred[i + 1] / 2 + pressure_pred[i + 1] / GlobalDefs.BLOOD_DENSITY + g_energy[i + 1]) - (velocity_pred[i] * velocity_pred[i] / 2 + pressure_pred[i] / GlobalDefs.BLOOD_DENSITY + g_energy[i]));
               velocity[i] = velocity[i] - 1.0 / 2.0 / GlobalDefs.BLOOD_DENSITY / lumen_sq[i] * (GlobalDefs.FRICTION_C * viscosity[i] * Math.PI * velocity[i]) * dt;
               lumen_sq[i] = (lumen_sq[i] + lumen_sq_pred[i]) / 2 - dt / dz / 2 * (lumen_sq_pred[i + 1] * velocity_pred[i + 1] - lumen_sq_pred[i] * velocity_pred[i]);
               pressure[i] = calcPressure(lumen_sq[i],   i);

           }
            //for (int i = 0; i < L; i++)
            //{
            //    if (nodes[i].blocked == true)
            //        velocity[i] = 0;
            //}

            if (L == 2)
           {
               velocity[0] = (velocity[0] + velocity[1])/2; 
               velocity[1] = velocity[0];

               pressure[0] = (pressure[0] + pressure[1]) / 2;
               pressure[1] = pressure[0];

               lumen_sq[0] = (lumen_sq[0] + lumen_sq[1]) / 2;
               lumen_sq[1] = lumen_sq[0];
           }

           if (L > 2)
           {
               
               velocity[L - 1] = velocity[L - 2];
               velocity[0] = velocity[1];

               pressure[L - 1] = pressure[L - 2];
               pressure[0] = pressure[1];
                               
               lumen_sq[L - 1] = lumen_sq[L - 2];
               lumen_sq[0] = lumen_sq[1];
           }        
#endif
            this.state2nodes();
        }

        protected double calcPressure(double _lumen_sq, int i)
        {
            return GlobalDefs.DIASTOLIC_PRESSURE + (beta_1[i] * (Math.Sqrt(_lumen_sq) - Math.Sqrt(lumen_sq_0[i])));
        } 
       
        protected double calcPressure(int i)
        {
            return GlobalDefs.DIASTOLIC_PRESSURE + (beta_1[i] * (Math.Sqrt(lumen_sq[i]) - Math.Sqrt(lumen_sq_0[i])));
        }
        
        protected double calcLumen_sq(int i)
        {
            return (double)Math.Pow((pressure[i] - GlobalDefs.DIASTOLIC_PRESSURE) / beta_1[i] + Math.Sqrt(lumen_sq_0[i]), 2);
        }


        protected double young_modulus;
        protected double[] wall_thickhess;
        protected double[]      viscosity;
        protected double   beta;
        protected double[] beta_1;
        protected double[] lumen_sq_0;
        protected double[] g_energy;


        protected double curr_dt;

#if FAST
        protected McCormackThread thread_solver;
#endif
     
    }

/*   public class ViscoElasticThread : ElasticThread
    {        
        public ViscoElasticThread(Thread _protothread, GetBetaFunction getElsticBeta): base(_protothread, getElsticBeta)
        {
            int L = nodes.GetLength(0);
            
            unsafe
            {
                double* dZ = (double*)Marshal.AllocHGlobal((L - 1) * sizeof(double));
                for (int i = 0; i < L - 1; i++)
                    dZ[i] = Vector3.Distance(nodes[i + 1].position, nodes[i].position);


                thread_solver = new McCormackThread(L, dZ, GlobalDefs.BLOOD_DENSITY, GlobalDefs.BLOOD_VISC);

                for (int i = 0; i < L; i++)
                {
                    int I = i;
                    _1DFunction_del pressure_func_del = delegate(double x)
                    {
                        return calcPressure(x, I);
                    };

                    delegate1DFunc f = new delegate1DFunc(pressure_func_del);

                    thread_solver.addFunc(f, i);
                }

                fixed (double* g_ptr = &g_energy[0])
                {
                    thread_solver.setGravity(g_ptr);
                }
            }

            clotes = new List<Clot>(); 
        }

        unsafe public override void calcThread(double dt)
        {
            current_time += dt;
            double dz = 0;
            curr_dt = dt;

            if(nodes.First().id==1742)
                nodes.First().id = 1742;

            if (nodes.Last().id == 1742)
                nodes.Last().id = 1742;

            int L = nodes.GetLength(0);

            this.nodes2state();

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

            if (L > 1)
            {

                velocity[L - 1] = velocity[L - 2];
                velocity[0] = velocity[1];

                pressure[L - 1] = pressure[L - 2];
                pressure[0] = pressure[1];

                lumen_sq[L - 1] = lumen_sq[L - 2];
                lumen_sq[0] = lumen_sq[1];
            }

            this.state2nodes();
        }

       protected double calcPressure(double _lumen_sq, int i)
        {
           double Gamma = 2.0 / 3.0 * Math.Sqrt(Math.PI) * wall_thickhess[i] * GlobalDefs.phi;
           return Gamma / lumen_sq_0[i] / Math.Sqrt(_lumen_sq) * (_lumen_sq - lumen_sq_old[i]) / curr_dt + base.calcPressure(_lumen_sq, i);          
        }

        protected double curr_dt;
    }*/
}
    