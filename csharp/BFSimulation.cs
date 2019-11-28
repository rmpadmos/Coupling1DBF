using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Text.RegularExpressions;
using System.Diagnostics;
//using System.IO.MemoryMappedFiles;

namespace BloodFlow
{

	public enum InletType { FLUX, PRESSURE }
	public enum SolutionState { CONTINUE, ERROR, FINISHED }

	class GlobalDefs
	{
        // overwrite these with values defined in the parameter file
		//public static float YOUNG_MODULUS = 225.0e3f;//Pa        
		public static double BLOOD_DENSITY = 1040f; //kg/m^3 
		public static double BLOOD_VISC = 3.5e-3f; //Pa*s
		public static double DIASTOLIC_PRESSURE = 10.0e+3f; //Pa, taken as reference and initial guess
		public static double OUT_PRESSURE = 0.0f;
        public static double FRICTION_C = 22.0f;
        public static string ScriptLocation = "./../../../../scripts/Check_Convergence.py";

        // not used at the moment
        //public static double DIASTOLIC_PRESSURE_1 = 0;
        //public static double SISTOLIC_PRESSURE = 15.9e+3f; //Pa
        public static double GRAVITY = 0.0f;// m/s
		public static double HEART_PERIOD = 1.0;
        public static double phi = 0;
		public static Vector1 ZERO_POINT = new Vector1(0.0f);
		public static Vector1 DOWN = new Vector1(0.0f);
		public static int av_iter_num = 0;

		static public TableFunction agent_inlet_flux = delegate (double a)
		{
			return 0.0;
		};

		static public double getBoileauBeta(double R0, double elasticity)
		{
			double h_a = 0.2802;
			double h_b = -505.3; //m^-1
			double h_c = 0.1324;
			double h_d = -11.14; //m^-1 
			double y_m = elasticity;

			double w_t = R0 * (h_a * Math.Exp(h_b * R0) + h_c * Math.Exp(h_d * R0));
			return 4.0 / 3.0 * Math.Sqrt(Math.PI) * y_m * w_t;
		}

		static public double getFixedHBeta(double R0, double elasticity)
		{
			double w_t = getFixedWallThickness(R0);
			return 4.0 / 3.0 * Math.Sqrt(Math.PI) * elasticity * w_t;
		}

		static public double getFixedWallThickness(double R0)
		{
			return 1.5e-3;
		}

		static public double getBoileauWallThickness(double R0, double elasticity)
		{
			double h_a = 0.2802;
			double h_b = -505.3; //m^-1
			double h_c = 0.1324;
			double h_d = -11.14; //m^-1 
			double y_m = elasticity;

			return R0 * (h_a * Math.Exp(h_b * R0) + h_c * Math.Exp(h_d * R0));
		}

		static public double getSoundSpeed(double R0, double elasticity)
		{
			double beta = getBoileauBeta(R0, elasticity);
			double A0 = Math.PI * R0 * R0;
			return Math.Sqrt(beta / (2.0 * BLOOD_DENSITY)) * Math.Pow(A0, -0.25);
		}

		static public double getR1(double R0, double elasticity)
		{
			double A0 = Math.PI * R0 * R0;
			return (BLOOD_DENSITY * getSoundSpeed(R0, elasticity)) / A0;
		}

		static public double sin_flux(double t)
		{
			return (double)Math.Sin(Math.PI * t * 10) * 0.25f + 10f; // ml/s
		}

		static public double single_pusle(double t)
		{
			return (double)1.0e-7 * Math.Exp(-1.0e4 * (t - 0.05) * (t - 0.05)); // ml/s
		}
	}

	public struct BodyPartInlet
	{
		public BodyPartInlet(string _name, List<VascularNode> _inlet, List<double> _inlet_flux)
		{
			name = _name;
			inlet = new List<VascularNode>(_inlet);
			inlet_flux = new List<double>(_inlet_flux);
			outlet = new List<BoundaryCondition>();
		}
		public string name;
		public List<VascularNode> inlet;
		public List<BoundaryCondition> outlet;
		public List<double> inlet_flux;
	}

	public class BFSimulator
	{
		public BFSimulator(string top_text, int[] inlet_nodes, List<ClotNode> clotnodes)
		{
			timestep_N = 0;
			init = true;

			end_time = -1.0f;
			current_time = 0.0f;

			IO_Module.localization = "en-US";

		//TableFunction heart_inlet = null;
			List<BC_params> rcr_params = new List<BC_params>();
			try
			{
				IO_Module.LoadTopologyFromString(top_text, out v_net);
			}
			catch { }


			GlobalDefs.ZERO_POINT = v_net.vascular_system[0].position;

			getFloatValueDelegate getProximaDst; // special delerates for wide-width search on arterial network
			setFloatValueDelegate setProximaDst; // -"-

			v_net.defineNodeDirVectors(inlet_nodes, out getProximaDst, out setProximaDst);

			v_net.defineNet(getProximaDst, setProximaDst, clotnodes);

			for (int i = 0; i < v_net.threads.Count; i++)
				v_net.specifyThreadType(i, new ElasticThread(v_net.threads[i], GlobalDefs.getBoileauBeta));

			// WaveTransmissive BCs are default, just pass pressure wave outside (have very small reflection coeff.)
			// No inlet BC are by def. It should be defined separatelly.
			for (int i = 0; i < v_net.bounds.Count; i++)
				// v_net.bounds[i] = new WaveTransmissive(v_net.bounds[i], GlobalDefs.getBoileauBeta);//new PressureOutletRCR(v_net.bounds[i], GlobalDefs.getBoileauBeta, 0.1e+9, 0.1e+9, 1.0e-12); // //
				v_net.bounds[i] = new PressureOutletRCR(v_net.bounds[i], GlobalDefs.getBoileauBeta, 2.0e+9, 100.0e+9, 0.1e-12); // //

			for (int i = 0; i < v_net.knots.Count; i++)
				v_net.knots[i] = new StandartKnot(v_net.knots[i], GlobalDefs.getBoileauBeta);


			control_point_list = new List<NodeSummary>();

			state = SolutionState.CONTINUE;
		}

		public void setTimeParameters(float timestep, float end_time, float av_pariod)
		{
			this.delta_tau = timestep;
			this.end_time = end_time;
			this.stenosis_av_period = av_pariod;
		}

		public BodyPartInlet createBodyPartInlet(string _name, int[] _inlet, double[] _fluxes)
		{
			List<VascularNode> _inlet_nodes = new List<VascularNode>();
			for (int i = 0; i < _inlet.Length; i++)
				_inlet_nodes.Add(v_net.vascular_system.Find(x => x.id == _inlet[i]));

			return new BodyPartInlet(_name, _inlet_nodes, _fluxes.ToList());

		}

		public BFSimulator(string path, string top_filename, int[] inlet_nodes, List<ClotNode> Clot_task)
			: this(
			File.ReadAllText(path + top_filename), inlet_nodes, Clot_task
           )
		{

		}

		public void setInletBC(string inlet_data, int node_number, InletType type)
		{
			TableFunction inlet_function = null;
			IO_Module.LoadTableFunctionFromString(inlet_data, out inlet_function);
			int id = v_net.bounds.FindIndex(x => x.core_node.id == node_number);

			if (type == InletType.FLUX)
				v_net.bounds[id] = new InletFlux(v_net.bounds[id], inlet_function, GlobalDefs.getBoileauBeta);

			if (type == InletType.PRESSURE)
				v_net.bounds[id] = new InletPressure(v_net.bounds[id], inlet_function, GlobalDefs.getBoileauBeta);
		}

		public void setOutletBC(string parameters)
		{
			List<PressureOutletRCR> oultlet_bc_set = new List<PressureOutletRCR>();

			if (parameters != "")
			{
				List<BC_params> rcr_params = new List<BC_params>();
				IO_Module.LoadBC_RCR_paramsFromString(parameters, out rcr_params);

				foreach (var bc_params in rcr_params)
				{
					try
					{
						int ind = v_net.bounds.FindIndex(x => x.core_node.id == bc_params.id);
						v_net.bounds[ind] = new PressureOutletRCR(v_net.bounds[ind], GlobalDefs.getBoileauBeta, bc_params.R1, bc_params.R2, bc_params.C);
						oultlet_bc_set.Add((PressureOutletRCR)v_net.bounds[ind]);
					}
					catch { }
				}

				return;
			}

			for (int i = 0; i < v_net.bounds.Count; i++)
			{
				BoundaryCondition bn = v_net.bounds[i];
				bn = new PressureOutletRCR(bn, GlobalDefs.getBoileauBeta, 0.0, 0.0, 0.0);
				oultlet_bc_set.Add((PressureOutletRCR)bn);
			}
		}

		public void setOutletBC(int start_point_id, List<BodyPartInlet> parts, float Q_total, float C_total, float R_total)
		{
			Dictionary<int, List<VascularNode>> outlet_groups = new Dictionary<int, List<VascularNode>>();

			List<VascularNode> reference_points = new List<VascularNode>();
			foreach (var p in parts)
				reference_points.AddRange(p.inlet);

			List<VascularNode> front_0 = new List<VascularNode>();
			List<VascularNode> front_1 = new List<VascularNode>();

			getFloatValueDelegate get_dist;
			setFloatValueDelegate set_dist;

			getFloatValueDelegate get_group;
			setFloatValueDelegate set_group;

			VascularNode.newFloatValueLayer(out get_dist, out set_dist);
			VascularNode.newFloatValueLayer(out get_group, out set_group);

			VascularNode start_point = v_net.vascular_system.Find(x => x.id == start_point_id);

			front_0.Add(start_point);
			set_dist(start_point, 0);

			/*  Q_total = 0;
			  for (int i = 0; i < parts.Count; i++)
				  Q_total += (float)parts[i].inlet_flux.Sum();*/

			while (true)
			{
				foreach (var n in front_0)
				{
					if (reference_points.Contains(n) && (!outlet_groups.ContainsKey(n.id)))
					{
						outlet_groups.Add(n.id, new List<VascularNode>());
						set_group(n, n.id);
					}

					foreach (var nn in n.neighbours)
						if (get_dist(nn) > get_dist(n))
						{
							set_dist(nn, get_dist(n) + 1);
							if (!front_1.Contains(nn))
							{
								front_1.Add(nn);
								set_group(nn, get_group(n));
								if (outlet_groups.ContainsKey((int)get_group(n)))
									outlet_groups[(int)get_group(n)].Add(nn);
							}
						}
				}

				front_1 = front_1.Distinct().ToList();
				front_0 = new List<VascularNode>(front_1);
				front_1.Clear();

				if (front_0.Count == 0)
					break;
			}

			foreach (var part in parts)
			{
				foreach (var root in part.inlet)
				{
					if (outlet_groups.ContainsKey(root.id))
						foreach (var n in outlet_groups[root.id])
							if (n.neighbours.Count == 1)
							{
								part.outlet.Add(v_net.bounds.Find(x => x.core_node == n));
							}
				}
			}

			for (int i = 0; i < parts.Count; i++)
			{
				double outlet_cube_summ = 0;
				double inlet_cube_summ = 0;
				foreach (var on in parts[i].outlet)
				{
					outlet_cube_summ += Math.Pow(on.core_node.radius, 3);
				}

				foreach (var inp_n in parts[i].inlet)
				{
					inlet_cube_summ += Math.Pow(inp_n.radius, 3);
				}

				Console.WriteLine(parts[i].name + "    inlet/outlet: " + inlet_cube_summ / outlet_cube_summ);
			}

			for (int i = 0; i < parts.Count; i++)
			{
				double a_cube_outlet = 0;
				double tot_inlet_flux = 0;

				for (int j = 0; j < parts[i].outlet.Count; j++)
					a_cube_outlet += Math.Pow(parts[i].outlet[j].core_node.radius, 3);


				for (int j = 0; j < parts[i].inlet_flux.Count; j++)
					tot_inlet_flux += parts[i].inlet_flux[j];

				double R_part = (Q_total / tot_inlet_flux) * R_total;

				foreach (PressureOutletRCR bc in parts[i].outlet)
				{
					double Rt = a_cube_outlet / Math.Pow(bc.core_node.radius, 3) * R_part;
					// inv_tot_R += 1.0 / Rt;
					double R2 = Rt - bc.R1;
					if (R2 < 0)
					{
						R2 = 0.1;
						bc.R1 = Rt - R2;
					}

					bc.R2 = R2;
					bc.C = C_total * R_total / Rt;
				}
			}
		}

		public void setOutletBC(List<BC_params> rcr_params)
		{
			List<PressureOutletRCR> oultlet_bc_set = new List<PressureOutletRCR>();

			foreach (var bc_params in rcr_params)
			{
				try
				{
					int ind = v_net.bounds.FindIndex(x => x.core_node.id == bc_params.id);
					v_net.bounds[ind] = new PressureOutletRCR(v_net.bounds[ind], GlobalDefs.getBoileauBeta, bc_params.R1, bc_params.R2, bc_params.C);
					oultlet_bc_set.Add((PressureOutletRCR)v_net.bounds[ind]);
				}
				catch { }
			}
		}

		public void Update()
		{
			Update(delta_tau);
		}

        public void ResetTime()
        {
            timestep_N = 0;
            current_time = 0;

            Parallel.ForEach(v_net.threads, (tr) =>
            {
                tr.current_time = 0;
            });
            Parallel.For(0, (v_net.bounds.Count - Program.Coupled_Nodes_N), index =>
            {
                v_net.bounds[index].current_time=0;
            });
            Parallel.ForEach(v_net.knots, (kn) =>
            {
                kn.current_time = 0;
            });
        }

        public void Update( float timestep)
		{
			this.delta_tau = timestep;
            timestep_N++;
            current_time = timestep_N * delta_tau;
			
			Parallel.ForEach(v_net.threads, (tr) =>
			{
				tr.calcThread(delta_tau);
			});//*/

            //  foreach( var tr in v_net.threads)
            //    tr.calcThread(delta_tau);

            //Parallel.ForEach(v_net.bounds, (bc) =>
            //{
            //	bc.doBC(delta_tau);
            //});                    //*/
            //foreach (var bc in v_net.bounds)
            //{
            //	bc.doBC(delta_tau);
            //}//*/



            Parallel.For(0, (v_net.bounds.Count - Program.Coupled_Nodes_N), index =>
			{
				v_net.bounds[index].doBC(delta_tau);//do not doBC on steadyBC nodes. (iterative function)
													//Console.WriteLine(index);
			});

            // collect data of the other BC nodes
            if (Program.Coupled_Nodes_N > 0)
            {
    			for (int i = v_net.bounds.Count - Program.Coupled_Nodes_N; i < v_net.bounds.Count; i++)
    			{
					Program.stream.Write((i-v_net.bounds.Count+Program.Coupled_Nodes_N)*8, v_net.bounds[i].core_node.velocity * v_net.bounds[i].core_node.lumen_sq * v_net.bounds[i].v_sign[0]*1e9);
    			    //Console.WriteLine("Ping");
    			}
			    //	Console.WriteLine(Program.stream.ReadDouble(0));			
				waiting = true;
				Program.Process.StandardInput.WriteLine("1");//send data to the steady_network_solver, upon responce, update BC nodes.
				while (waiting == true)
				{
					System.Threading.Thread.Sleep(1);
					//Program.Process.WaitForInputIdle();
				}
			}

            //foreach (var kn in v_net.knots)
            //  kn.doCoupling(delta_tau);
            Parallel.ForEach(v_net.knots, (kn) =>
            {
                kn.doCoupling(delta_tau);
            });

            //Parallel.ForEach(v_net.threads, (tr) =>
            //{
            //	tr.updateState();
            //	// tr.updateStateFFR();
            //});

            //foreach (var clotnode in ClotNodes)
            //{
            //clotnode.Node.velocity = 0;
            //clotnode.Node.lumen_sq = clotnode.Node.lumen_sq_0;
            //}

            //foreach (var tr in v_net.threads)
            //{
            //	tr.updateState();
            //	tr.updateStateFFR();
            //}

            //foreach (var cp in control_point_list)
            //{
            //	cp.addFluxVal();
            //	cp.addPrsVal();
            //}
        }
		public void Set_Pressure_Steady_Solver(object sendingProcess,
		  DataReceivedEventArgs outLine)
		{
			//if (!String.IsNullOrEmpty(outLine.Data))
			//{
				//Console.WriteLine(outLine.Data);
				//Steady_Output = outLine.Data;
				//var Pressure_steady = outLine.Data.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
			Parallel.For(0, Program.Coupled_Nodes_N, index =>
				{
				v_net.bounds[v_net.bounds.Count - Program.Coupled_Nodes_N + index].doBCsteady(delta_tau, Program.stream.ReadDouble((index + Program.Coupled_Nodes_N)*8));
				});
				waiting = false;
			//}
		}
		//public void Set_Pressure_Steady_Solver(object sendingProcess,
    //      DataReceivedEventArgs outLine)
    //    {
			
    //        Parallel.For(0, 4, index =>
    //            {
				//v_net.bounds[v_net.bounds.Count - 4 + index].doBCsteady(delta_tau, Program.stream.ReadDouble(4+index));
        //        });
        //        waiting = false;            
        //}
		public SolutionState Control()
		{
			//solution_state = SolutionState.CONTINUE;

			//if (end_time > 0)
			if (current_time >= end_time)
					solution_state = SolutionState.FINISHED;

			if (current_time % stenosis_av_period < delta_tau)
				foreach (var cp in control_point_list)
					cp.reset();

			foreach (var bc in v_net.bounds)
			{
				if (bc.ValueControl() >= 0)
					state = SolutionState.ERROR;
			}

			foreach (var kn in v_net.knots)
			{
				if (kn.NaNControl() >= 0)
					state = SolutionState.ERROR;
			}
			solution_state = state;
			return state;
		}

		public void addContorlPoint(int node_id)
		{
			NodeSummary nd_summary = new NodeSummary();
			nd_summary.node = v_net.vascular_system.Find(x => x.id == node_id);
			control_point_list.Add(nd_summary);
		}

		public NodeSummary getContorlPoint(int node_id)
		{
			try
			{
				return control_point_list.Find(x => x.node.id == node_id);
			}
			catch
			{
				return null;
			}
		}

		public void FullReset()
		{
			v_net.fullReset();
			clotes_id.Clear();
			current_time = 0;

			foreach (var cp in control_point_list)
				cp.reset();

			solution_state = SolutionState.CONTINUE;
		}

		public bool removeLastClot()
		{
			if (!clotes_id.Any())
				return false;

			v_net.removeClot(clotes_id.Pop());
			return true;
		}

		public List<int> getClotesId()
		{
			return clotes_id.ToList();
		}

		public void removeAllClots()
		{
			while (removeLastClot())
			{ }
		}

		public void addClot(int node_id, float degree)
		{
			if (degree > 0.99)
				degree = 0.98f;

			if (v_net.setCloth(node_id, degree))
				clotes_id.Push(node_id);
		}

		public bool isClotSet()
		{
			if (clotes_id.Count == 0)
				return false;

			return true;
		}

		public IWriteVascularNet getVNet()
		{
			return v_net;
		}


		public List<int> getCenterNodesID(int min_len, float max_R)
		{
			List<int> center_nodes_id = new List<int>();
			foreach (var tr in this.v_net.threads)
				if (tr.nodes.GetLength(0) >= min_len && tr.nodes[tr.nodes.GetLength(0) / 2].radius < max_R * 1e-3)
					center_nodes_id.Add(tr.nodes[tr.nodes.GetLength(0) / 2].id);
			return center_nodes_id;
		}


		public List<NodeSummary> control_point_list;
		public SolutionState solution_state;
		protected VascularNet v_net;

		//  public NodeSummary[] thread_summary;

		protected bool init;
		public float delta_tau { get; protected set; }
		public int timestep_N { get; protected set; }
		public double current_time { get; protected set; }
		public float stenosis_av_period { get; protected set; }
		//       public float stabilisation_time { get; protected set; }       

		public float output_step { get; set; }
		public float end_time { get; set; }
		public bool waiting = true;
		protected Stack<int> clotes_id = new Stack<int>();

        public List<ClotNode> ClotNodes;

        SolutionState state;
		//private object process;
	}

}
