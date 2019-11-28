using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using System.Diagnostics;

namespace BloodFlow
{
    public struct BC_params
    {
        public int id;
        public double R1, R2, C;

        public BC_params(int id, double R1, double R2, double C)
        {
            this.id = id;
            this.R1 = R1;
            this.R2 = R2;
            this.C = C;
        }
    }

    public class BoundaryCondition
    {
        protected static float        TIMEGAP = 1e-3f;
        protected static float DEF_BASEPERIOD =  1.0f;

        public BoundaryCondition(VascularNode _core_node, double start_time)
        {
            core_node      = _core_node;
            //Only 1 neughbour is possible for terminal node
            neighbour_node = _core_node.neighbours[0];

            current_time  = start_time;
            previous_time = current_time - TIMEGAP;
            v_sign = new int[2];            
            DefineSign();
        }

        public virtual void reset()
        {
            current_time = 0;
            previous_time = current_time - TIMEGAP;
            core_node.velocity = 0;
            core_node.pressure = 0;
            core_node.lumen_sq = core_node.lumen_sq_0;
        }

        public virtual void doBC(double dt)
        {
            previous_time = current_time;
            current_time += dt;            
        }

		public virtual void doBCsteady(double dt, double pressure)
        {
            previous_time = current_time;
            current_time += dt;
        }


        public virtual double getTime()
        { return current_time; }

        public virtual int ValueControl()
        {
            if (double.IsNaN(core_node.pressure))
                return core_node.id;
            if (double.IsNaN(core_node.velocity))
                return core_node.id;
            return -1;
        }

        protected virtual void DefineSign()
        {
            Vector1 dir_vector1 = new Vector1();
            Vector1 dir_vector2 = new Vector1();

            // Positive direction is outflow from termianl node
            dir_vector1 = neighbour_node.position - core_node.position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));

            dir_vector2 = neighbour_node.dir_vector;
            v_sign[1] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));
        }

        public double current_time;
        public double previous_time;
        public double base_period;

        public int[] v_sign;
        public VascularNode      core_node;
        public VascularNode neighbour_node;
    };


    public class InletFlux : BoundaryCondition
    {
        public InletFlux(BoundaryCondition BC, TableFunction _flux, GetBetaFunction getElasticBeta)
            : base(BC.core_node, BC.current_time)
        {            
            this.previous_velocity = 0;         

            base_flux_on_time = _flux;
                 flux_on_time = _flux;


            // Artery lumen of terminal and previous nodes are set the same
            core_node.lumen_sq_0 = neighbour_node.lumen_sq_0;
            core_node.lumen_sq   = core_node.lumen_sq_0;                
            

            double R0 = Math.Sqrt(core_node.lumen_sq_0 / Math.PI);
            
            
            beta_1 = getElasticBeta(R0, core_node.elasticity) / core_node.lumen_sq_0;
            flux_function = delegate(double A)
            {
                return A * chrt_back + 4 * Math.Pow(A, 1.25f) * Math.Sqrt(beta_1 / 2.0 /  GlobalDefs.BLOOD_DENSITY) - flux_on_time(current_time);
            };

            this.chrt_back = -4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);            
        }

        public override void doBC(double dt)
        {
            if (true)
            {
                previous_time = current_time;
                current_time = current_time + dt;
                previous_velocity = current_velocity;

                double flux = flux_on_time(current_time);
                //Console.WriteLine(flux);

                double inlet_lumen = myMath.NewtonSolver(flux_function, core_node.lumen_sq, core_node.lumen_sq_0*1e-4f, core_node.lumen_sq_0 * 1e-6);
                double U = flux / inlet_lumen;


                double dZ = Vector1.Distance(core_node.neighbours.Last().position, core_node.position);
                double dUdt = (U - previous_velocity) / dt;
                double dUdZ = (core_node.neighbours.Last().velocity * v_sign[1] - U) / dZ;
                double visc_term = (double)(GlobalDefs.FRICTION_C * GlobalDefs.BLOOD_VISC * Math.PI * U / GlobalDefs.BLOOD_DENSITY / inlet_lumen);
                double P = core_node.neighbours.Last().pressure + (dUdt + U * dUdZ + visc_term) * dZ *  GlobalDefs.BLOOD_DENSITY;


                core_node.velocity = U * v_sign[0];
                core_node.pressure = P;
                core_node.lumen_sq = inlet_lumen;
                //Console.WriteLine(core_node.velocity);
                chrt_back = core_node.velocity - 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f /  GlobalDefs.BLOOD_DENSITY);                
                current_velocity = U;
            }
        }

        public override void reset()
        {
            base.reset();
            
            chrt_back = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);            
            current_velocity = 0;
            previous_velocity = 0;            
        }

        protected SimpleFunction flux_function;
        protected double beta_1;
        

        protected double previous_velocity;
        protected double current_velocity;
        protected double chrt_back;

        // Need for change heart rate////////////        
        public TableFunction flux_on_time;
        public readonly TableFunction base_flux_on_time;
        //public readonly float base_period;
        /////////////////////////////////////////

        // Not used /////////////
        protected double Q_min;
        protected double Q_max;
        protected double diastolic_pressure;
        protected double sistolic_pressure;

        protected Queue<double> pressure_hist;

        protected double pulse_time_interval;
        protected double flux_minmax_time_interval;

        //private float sample_dt;
        //////////////////////////
    };

    public class PressureOutletRCR : BoundaryCondition
    {
        public PressureOutletRCR(BoundaryCondition BC, GetBetaFunction getElsticBeta, double _R1, double _R2, double _C)
            : base(BC.core_node, BC.current_time)
        {
            
            // Windkessel params, R2 is parallel to C //
            R1 = _R1;
            R2 = _R2;
            C = _C;
            ////////////////////////////////////////////
            
            //tx_ij - i - time moment, j - space point, 1 - current moment (core node), 0 -previous moment (previus point)/
            U_tx_10 = 0;
            U_tx_01 = 0;
            P_tx_10 = GlobalDefs.DIASTOLIC_PRESSURE;
            P_tx_01 = GlobalDefs.DIASTOLIC_PRESSURE;
            P_tx_11 = GlobalDefs.DIASTOLIC_PRESSURE;
            A_tx_01 = core_node.lumen_sq_0;
            A_tx_10 = core_node.lumen_sq_0;
            Q_t_1 = 0;
            Q_t_0 = 0;
            ///////////////////////////////////////////

            beta_1 = GlobalDefs.getBoileauBeta(core_node.radius,core_node.elasticity) / core_node.lumen_sq_0;
            chrt_function = delegate(double A_tx_11)
            {
                double chrt_frw_right = Q_t_1 / A_tx_11 + 4 * Math.Pow(A_tx_11, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
                return chrt_frw_left - chrt_frw_right;
            };

            DefineSign();
        }

        protected override void DefineSign()
        {
            Vector1 dir_vector1 = new Vector1();
            Vector1 dir_vector2 = new Vector1();


            dir_vector1 = core_node.position - core_node.neighbours.First().position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));

            dir_vector2 = core_node.neighbours.Last().dir_vector;
            v_sign[1] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));
        }

        public override int ValueControl()
        {
            if (double.IsNaN(core_node.pressure))
                return core_node.id;
            return -1;
        }
        
        public override void doBC(double dt)
        {
             
            interfunc_dx = Vector1.Distance(core_node.position, core_node.neighbours.Last().position);
            interfunc_dt = dt;
            previous_time = current_time;
            current_time = current_time + dt;
			//Console.WriteLine(core_node.position);
            //if (core_node.id == 13234)
                //core_node.id = 13234;

			Q_t_1 = core_node.lumen_sq * core_node.velocity * v_sign[0];
            //tx_ij - i - time moment, j - space point, 1 - current moment (core node), 0 -previous moment (previus point)/
            // Euler scheme for ODE integration///
            if (core_node.nodetype >= 3)
            {
                Q_t_1 = 0.0;
            }
            double dQdt = (Q_t_1 - Q_t_0) / dt;//Console.WriteLine(dQdt);

			////////////////////////////////////
			// Send flow rate to 1-D steady model
			// Wait for answer
			// Set pressure 
            //string flags = " -Qin " + Q_t_1.ToString();
            //Process Process = new Process();

            //Process.StartInfo.UseShellExecute = false;
            //Process.StartInfo.FileName = "./steady_network_solver";
            //Process.StartInfo.Arguments = flags;
            //Process.StartInfo.CreateNoWindow = true;
            //Process.StartInfo.RedirectStandardOutput = true;
            //Process.Start();
            ////q = Process.StandardOutput.ReadLine();              
            //Process.WaitForExit();
            //string Pressure_reading = System.IO.File.ReadAllText("p_in.dat");
            //P_tx_11 = Convert.ToDouble(Pressure_reading);    
            
            P_tx_10 = neighbour_node.pressure;
            U_tx_10 = neighbour_node.velocity * v_sign[1];
            A_tx_10 = neighbour_node.lumen_sq;
            if (core_node.nodetype >= 3)
            {
                P_tx_11 = P_tx_10;
            }
            else
            {
                P_tx_11 = (Q_t_0 * (1 + R1 / R2) + C * R1 * dQdt + GlobalDefs.OUT_PRESSURE / R2 + P_tx_01 * C / dt) / (C / dt + 1 / R2);
            }
            //P_tx_11 = (Q_t_0 * (1 + R1 / R2) + C * R1 * dQdt + GlobalDefs.OUT_PRESSURE / R2 + P_tx_01 * C / dt) / (C / dt + 1 / R2);
            chrt_frw_left = U_tx_01 + 4 * Math.Pow(A_tx_10, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
            //if (core_node.nodetype >= 3)
            //{
            //    A_tx_11 = core_node.lumen_sq_0;
            //}
            //else
            //{
            //    A_tx_11 = myMath.NewtonSolver(chrt_function, A_tx_10, 1e-9, 1e-10);
            //}
            A_tx_11 = myMath.NewtonSolver(chrt_function, A_tx_10, 1e-9, 1e-10);
            U_tx_11 = Q_t_1 / A_tx_11;
                                   
            core_node.velocity = U_tx_11 * v_sign[0];
            core_node.lumen_sq = A_tx_11;
            core_node.pressure = P_tx_11;
            /*  core_node.neighbours.Last().pressure = core_node.pressure;
            core_node.neighbours.Last().velocity = core_node.velocity;  
            core_node.neighbours.Last().lumen_sq = core_node.lumen_sq;*/
            
            P_tx_01 = P_tx_11;
            A_tx_01 = A_tx_11;
            U_tx_01 = U_tx_11;
            Q_t_0 = Q_t_1;

        }

		public override void doBCsteady(double dt, double pressure)
        {

            interfunc_dx = Vector1.Distance(core_node.position, core_node.neighbours.Last().position);
            interfunc_dt = dt;
            previous_time = current_time;
            current_time = current_time + dt;
			//Console.WriteLine(R1);    
			//Console.WriteLine(R2);    
			//Console.WriteLine(C);  
            
			Q_t_1 = core_node.lumen_sq * core_node.velocity * v_sign[0];
            //tx_ij - i - time moment, j - space point, 1 - current moment (core node), 0 -previous moment (previus point)/
            // Euler scheme for ODE integration///
            //double dQdt = (Q_t_1 - Q_t_0) / dt;
            P_tx_11 = pressure;

            P_tx_10 = neighbour_node.pressure;
            U_tx_10 = neighbour_node.velocity * v_sign[1];
            A_tx_10 = neighbour_node.lumen_sq;

            chrt_frw_left = U_tx_01 + 4 * Math.Pow(A_tx_10, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);

            A_tx_11 = myMath.NewtonSolver(chrt_function, A_tx_10, 1e-9, 1e-10);
            U_tx_11 = Q_t_1 / A_tx_11;

            core_node.velocity = U_tx_11 * v_sign[0];
            core_node.lumen_sq = A_tx_11;
            core_node.pressure = P_tx_11;
            /*  core_node.neighbours.Last().pressure = core_node.pressure;
              core_node.neighbours.Last().velocity = core_node.velocity;
              core_node.neighbours.Last().lumen_sq = core_node.lumen_sq;*/

            P_tx_01 = P_tx_11;
            A_tx_01 = A_tx_11;
            U_tx_01 = U_tx_11;
            Q_t_0 = Q_t_1;
        }

        public double calcLumenSq(double pressure)
        {
            return (double)Math.Pow((pressure - GlobalDefs.DIASTOLIC_PRESSURE) / beta_1 + Math.Sqrt(core_node.lumen_sq_0), 2);
        }

        public double calcPressure(double lumen)
        {
            return GlobalDefs.DIASTOLIC_PRESSURE + beta_1 * (Math.Sqrt(core_node.lumen_sq) - Math.Sqrt(core_node.lumen_sq_0));
        }

        protected SimpleFunction chrt_function;

        protected double U_tx_01, U_tx_10, U_tx_11, P_tx_10, P_tx_11, P_tx_01;
        protected double A_tx_10, A_tx_01, A_tx_11, Q_t_0, Q_t_1, chrt_frw_left;
        protected double interfunc_dx, interfunc_dt;


        public double R1, R2, C, beta_1;

    }

    public class InletPressure : BoundaryCondition
    {
        public InletPressure(BoundaryCondition BC, TableFunction _pressure, GetBetaFunction getElasticBeta)
            : base(BC.core_node, BC.current_time)
        {            
            pressure_on_time = _pressure;

            double R0 = Math.Sqrt(core_node.lumen_sq_0 / Math.PI);
            beta_1 = getElasticBeta(R0, core_node.elasticity) / core_node.lumen_sq_0;

            
            core_node.lumen_sq_0 = neighbour_node.lumen_sq_0;
            core_node.lumen_sq = neighbour_node.lumen_sq_0;

            chrt = 4 * Math.Pow(core_node.lumen_sq_0, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);

            
        }

        public override void doBC(double dt)
        {
            previous_time = current_time;
            current_time = current_time + dt;

            double pressure = pressure_on_time(current_time);
            double inlet_lumen = Math.Pow((pressure - GlobalDefs.DIASTOLIC_PRESSURE) / beta_1 + Math.Sqrt(core_node.lumen_sq_0), 2);
            double neighbour_chrt = neighbour_node.velocity - 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25) * Math.Sqrt(beta_1 / GlobalDefs.BLOOD_DENSITY / 2.0f);
            double U = neighbour_chrt + 4 * Math.Pow(inlet_lumen, 0.25) * Math.Sqrt(beta_1 / GlobalDefs.BLOOD_DENSITY / 2.0f);

            core_node.velocity = U * v_sign[0];
            core_node.pressure = pressure;
            core_node.lumen_sq = inlet_lumen;
            chrt = core_node.velocity + 4 * Math.Pow(core_node.neighbours.Last().lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
            
            current_velocity = U;
        }

        protected TableFunction pressure_on_time;

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

        //private float sample_dt;
    };

    public class WaveTransmissive : BoundaryCondition
    {
        public WaveTransmissive(BoundaryCondition BC, GetBetaFunction getElsticBeta)
            : base(BC.core_node, BC.current_time)
        {           

            //tx_ij - i - time moment, j - space point, 1 - current moment (core node), 0 -previous moment (previus point)/
            Q_tx_10 = 0;
            Q_tx_01 = 0;
            P_tx_10 = GlobalDefs.DIASTOLIC_PRESSURE;
            P_tx_01 = GlobalDefs.DIASTOLIC_PRESSURE;
            P_tx_11 = GlobalDefs.DIASTOLIC_PRESSURE;            
            Q_t_1 = 0;
            Q_t_1 = 0;
            ///////////////////////////////////////////

            beta_1 = GlobalDefs.getBoileauBeta(core_node.radius, core_node.elasticity) / core_node.lumen_sq_0;

            chrt_function = delegate(double A_tx_11)
            {
                double chrt_frw_right = Q_t_1 / A_tx_11 + 4 * Math.Pow(A_tx_11, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
                return chrt_frw_left - chrt_frw_right;
            };


            DefineSign();
        }

        protected override void DefineSign()
        {
            Vector1 dir_vector1 = new Vector1();
            Vector1 dir_vector2 = new Vector1();


            dir_vector1 = core_node.position - core_node.neighbours.First().position;
            dir_vector2 = core_node.dir_vector;

            v_sign[0] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));

            dir_vector2 = core_node.neighbours.Last().dir_vector;
            v_sign[1] = Math.Sign(Vector1.Dot(dir_vector1, dir_vector2));
        }

        public override int ValueControl()
        {
            if (double.IsNaN(core_node.pressure))
                return core_node.id;
            return -1;
        }

        public override void doBC(double dt)
        {

            interfunc_dx  = Vector1.Distance(core_node.position, core_node.neighbours.Last().position);
            interfunc_dt  = dt;
            previous_time = current_time;
            current_time  = current_time + dt;


            Q_tx_01 = core_node.lumen_sq      * core_node.velocity      * v_sign[0];
            Q_tx_10 = neighbour_node.lumen_sq * neighbour_node.velocity * v_sign[1];
            U_tx_10 = neighbour_node.velocity * v_sign[1];

            P_tx_01 = core_node.pressure;
            P_tx_10 = neighbour_node.pressure;

            double sound_speed = Math.Pow(core_node.lumen_sq, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
            
            P_tx_11 = (interfunc_dx / sound_speed / dt * P_tx_01 - P_tx_10) / (1 + interfunc_dx / dt / sound_speed);
            Q_tx_11 = (interfunc_dx / sound_speed / dt * Q_tx_01 - Q_tx_10) / (1 + interfunc_dx / dt / sound_speed);                

            A_tx_10 = neighbour_node.lumen_sq;
            chrt_frw_left = U_tx_10 + 4 * Math.Pow(A_tx_10, 0.25f) * Math.Sqrt(beta_1 / 2.0f / GlobalDefs.BLOOD_DENSITY);
          //  A_tx_11 = myMath.NewtonSolver(chrt_function, A_tx_10, 1e-9, 1e-10);
            A_tx_11 = calcLumenSq(P_tx_11);


            core_node.velocity = Q_tx_11 / A_tx_11 * v_sign[0];
            core_node.lumen_sq = A_tx_11;
            core_node.pressure = P_tx_11;

            //A_tx_11 = calcLumenSq(P_tx_11);            
            
            Q_t_0 = Q_t_1;

        }

        public double calcLumenSq(double pressure)
        {
            return (double)Math.Pow((pressure - GlobalDefs.DIASTOLIC_PRESSURE) / beta_1 + Math.Sqrt(core_node.lumen_sq_0), 2);
        }

        public double calcPressure(double lumen)
        {
            return GlobalDefs.DIASTOLIC_PRESSURE + beta_1 * (Math.Sqrt(core_node.lumen_sq) - Math.Sqrt(core_node.lumen_sq_0));
        }

        protected SimpleFunction chrt_function;

        protected double U_tx_10, Q_tx_01, Q_tx_10, Q_tx_11, P_tx_10, P_tx_11, P_tx_01;
        protected double A_tx_11, A_tx_10,  Q_t_0, Q_t_1, chrt_frw_left, chrt_frw_right;
        protected double interfunc_dx, interfunc_dt;

        public double R1, R2, C, beta_1;

    }
}
