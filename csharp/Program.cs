using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using System.IO.MemoryMappedFiles;

namespace BloodFlow
{
    public enum adjustMode {AbsInit, AllInit, All, Personal, None};
    public enum CalculationMode {None, Stabilisation, WriteState, ReferenceAveraging, GetReferenceFlux, ClotAveraging, GetClotFlux, ClotStabilisations, RemoveClot, AddClot, ResetAveraging, FullResetAvergaing};    
    
    class Program
    {
		public static float  TIMESTEP =   0.1e-4f;//time between iterations
        public static float   AV_TIME =     1.0f;
        public static float  END_TIME =    12.0f;// simulate for 12 seconds
        //public static float WRITE_TIME =    10.0f;
		public static int SimulationSteps = (int)(END_TIME / TIMESTEP) + 1;//total iterations in the sim

        // stuff for communication to the perfusion model
		public static MemoryMappedViewAccessor stream;
		public static Process Process= new Process();
		public static int Coupled_Nodes_N;

		static void Main(string[] args)
		{
            // start timer
			Stopwatch stopWatch = new Stopwatch();
			stopWatch.Start();

			//BFSimulator bf_simuation = new BFSimulator("", @"verefication\full_body_Boileau_LR_2.5mm.top", GlobalDefs.getBoileauBeta);
			IO_Module.localization = "en-US";
			string top_filename = "";
			string parameters_filename = "";
			string out_filename = "";
			string stream_out_filename = "";
            string clot_filename = "";

            List<ClotDscr> stenosis_task = new List<ClotDscr>();
			List<Tuple<int, string>> inlet_data = new List<Tuple<int, string>>();

			string base_path = IO_Module.readTaskFile(args[0], ref top_filename, ref parameters_filename, ref inlet_data, ref stenosis_task, ref out_filename, ref Coupled_Nodes_N);

			string dyn_out_filename = @"Results.dyn";

            if (args.GetLength(0) > 1)
			{
				//out_filename = args[1] + out_filename;
				//dyn_out_filename = args[1] + dyn_out_filename;
                dyn_out_filename = args[1];
            }

            string ResultFolder = new FileInfo(dyn_out_filename).Directory.FullName+ Path.DirectorySeparatorChar;
            string ConvFile = ResultFolder + "Conv.txt";
            //List<VascularNode> ClotNodes = new List<VascularNode>();
            List<ClotNode> Clot_task = new List<ClotNode>();
            if (args.GetLength(0) > 2)
            {
                clot_filename = args[2];
                IO_Module.ReadClotFile(clot_filename, ref Clot_task);
                //foreach (var n in Clot_task):
                    //n.1

                //for (int i = 0; i < Clot_task.Count; i++)
                    //ClotNodes = bf_simuation.getVNet().vascular_system.Find(x => x.id == Clot_task[i].Node_id);
            }

            System.IO.StreamWriter out_dynamics_file = new System.IO.StreamWriter(dyn_out_filename);
			stream_out_filename = out_filename + "s";

			bf_simuation = new BFSimulator("", base_path + top_filename, inlet_data.Select(t => t.Item1).ToArray(), Clot_task);
            
            foreach (var inlt in inlet_data)
				bf_simuation.setInletBC(File.ReadAllText(base_path + inlt.Item2), inlt.Item1, InletType.FLUX);
           

            if (parameters_filename != "")//outlet parameters
				bf_simuation.setOutletBC(File.ReadAllText(base_path + parameters_filename));

            bf_simuation.setTimeParameters(TIMESTEP, END_TIME, AV_TIME);
            // Start the Steady_network_solver
            //Process Process = new Process();
            if (Coupled_Nodes_N > 0)
			{
				var mmf = MemoryMappedFile.CreateFromFile("/tmp/sharedfile", FileMode.OpenOrCreate, "/tmp/sharedfile", Coupled_Nodes_N * 16);
				stream = mmf.CreateViewAccessor();
				Process.StartInfo.UseShellExecute = false;
				//Process.StartInfo.FileName = "./steady_network_solver_shared";
				Process.StartInfo.FileName = "./coupled_solver";
				Process.StartInfo.Arguments = "-term_comm 0";
				Process.StartInfo.CreateNoWindow = true;
				Process.StartInfo.RedirectStandardOutput = true;
				Process.StartInfo.RedirectStandardInput = true;
				Process.Start();
				Console.WriteLine("Starting the Coupled Solver.");
				System.Threading.Thread.Sleep(20000);
				Process.BeginOutputReadLine();
				Process.OutputDataReceived += bf_simuation.Set_Pressure_Steady_Solver;
			}
            Process ConvCheck = new Process();
            ConvCheck.StartInfo.FileName = GlobalDefs.ScriptLocation;
            ConvCheck.StartInfo.Arguments = ResultFolder;
            ConvCheck.StartInfo.CreateNoWindow = true;

            //for (int i = 0; i < Clot_task.Count; i++)
            //{
            //    Clot_task[i].Node = bf_simuation.getVNet().vascular_system.Find(x => x.id == Clot_task[i].Node_id);
            //    Clot_task[i].Node.blocked = true;
            //}
            //bf_simuation.ClotNodes = Clot_task;

            int HBLoop = 0;
            while (true)
            {
                Console.WriteLine("Current Heartbeat loop: " + HBLoop);
                while (true)
                {
                    // Output and writing timers
                    if (CurrentIter % OUTPUT_PERIOD == 0 || CurrentIter == SimulationSteps)
                    {
                        Console.WriteLine("Time: " + bf_simuation.current_time + " Iteration: " + CurrentIter);
                        IO_Module.WriteState(bf_simuation.current_time, bf_simuation.getVNet(), out_dynamics_file);
                        //Console.WriteLine("State is written.");
                        //output is given as flowrate in ml/s, pressure in pa, and radius in mm
                        if ((bf_simuation.current_time % OUTPUT_PERIOD * 10) < bf_simuation.delta_tau)
                        {
                            out_dynamics_file.Flush();
                        }
                    }
                    if (bf_simuation.solution_state == SolutionState.ERROR)
                    {
                        Console.WriteLine("Error; physical time: " + bf_simuation.current_time);
                        break;
                    }

                    if (bf_simuation.solution_state == SolutionState.FINISHED)
                    {
                        Console.WriteLine("Physical time: " + bf_simuation.current_time);
                        Console.WriteLine("End");
                        break;
                    }
                    if (CurrentIter == Program.SimulationSteps)
                    {
                        break;
                    }

                    CurrentIter += 1;
                    bf_simuation.Control();// check for errors
                    bf_simuation.Update();// compute next timestep
                    //clot_timer += bf_simuation.delta_tau;
                }

                //check convergence
                out_dynamics_file.Close();
                ConvCheck.Start();               
                ConvCheck.WaitForExit();
                int Conv = 0;
                string[] lines = System.IO.File.ReadAllLines(ConvFile);
                Conv = int.Parse(lines[0]);      
                // if converged, stop
                if (Conv == 1 || HBLoop==20)
                    break;               
                out_dynamics_file = new System.IO.StreamWriter(dyn_out_filename);
                HBLoop += 1;
                // initiation of new heatbeat
                CurrentIter = 0;
                bf_simuation.ResetTime();
                //Console.WriteLine("Time: " + bf_simuation.current_time + " Iteration: " + CurrentIter);
                //IO_Module.WriteState(bf_simuation.current_time, bf_simuation.getVNet(), out_dynamics_file);
                //bf_simuation.timestep_N = 0;
                //bf_simuation.current_time = 0;
            }
            
			if (Coupled_Nodes_N > 0)// close program if it was open
			{
			Process.StandardInput.WriteLine("0");//close program
			Process.WaitForExit();
		    }
            // Output Elapsed time
            stopWatch.Stop();
            Console.WriteLine("Time elapsed: {0}", stopWatch.Elapsed);
        }
		public static int CurrentIter = 0;
		public static float clot_set_time = 0.0f;
        public static float STABILISATION_TIME = 10.0f;//10
        public static float AVERGAGE_PERIOD    = 2.0f;
        public static float CLOT_RELAXATION_PERIOD = 10.0f;//10
        public static float CLOT_REMOVE_RELAXATION_PERIOD = 1.0f;
        public static float OUTPUT_PERIOD = 0.02f;
        public static BFSimulator bf_simuation;
    }
}
