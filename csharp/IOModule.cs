using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Text.RegularExpressions;

namespace BloodFlow
{
    public delegate double TableFunction(double arg);

    class TopReadingException: System.Exception
    {
        public TopReadingException(string massage)
            : base(massage)
        { }
    }

    public class NodeSummary
    {
        public void reset()
        {
            summ_abs_flux = 0;
            num_of_mesurements_p = 0;
            num_of_mesurements_f = 0;
            Pd = double.MaxValue;
            Ps = 0;
            P_av = 0;           
        }

        public double get_av_flux()
        {
            if (num_of_mesurements_f != 0)
                return summ_abs_flux / num_of_mesurements_f;
            else return 0;
        }

        public double get_av_prs()
        {
            if (num_of_mesurements_p != 0)
                return P_av / num_of_mesurements_p;
            else return 0;
        }

        public void addFluxVal(float flx)
        {            
            latestFlux = flx;
            summ_abs_flux += flx;
            num_of_mesurements_f++;
        }

        public void addFluxVal()
        {
            addFluxVal((float)(node.velocity * node.lumen_sq));
        }

        public void addPrsVal(float prs)
        {
            P_av += prs;
            num_of_mesurements_p++;

            if (prs > Ps)
                Ps = prs;

            if (prs < Pd)
                Pd = prs;
        }

        public void addPrsVal()
        {
            addPrsVal((float)node.pressure);            
        }        

      public   VascularNode node;      
      public    int num_of_mesurements_p;
      public    int num_of_mesurements_f;
      public double proximal_dst;
      public double Pd;
      public double Ps;
      public double latestFlux;

      private double summ_abs_flux;
      private double P_av;

    }
    public class ClotNode
    {
        public int Clot_id;
        public VascularNode Node; 
        public int Node_id;
        public float Permeability;
        public float Porosity;
    }

        public class ClotDscr
    {
        public int   node_id;
        public float degree;
        public float reference_flux;
        public float depressed_flux;
        public bool is_ready = false;

        public float getRelDcr()
        {
            return depressed_flux / reference_flux;
        }
    }

    public class IO_Module
    {
        static public bool LoadTableFunctionFromFile(string filename, out TableFunction t_f)
        {
            return LoadTableFunctionFromString(File.ReadAllText(filename), out t_f);
        }


        static public bool LoadTopologyFromFile(string filename, out VascularNet r_vnet)
        {
            try
            {
                var text = File.ReadAllText(filename);
				return LoadTopologyFromString( text, out r_vnet );
            }
            catch
            {
                throw new TopReadingException("The topology file " + filename + " doesn't exist or can't be opened!\n" );              
            }
		}        

        static public bool LoadTopologyFromString(string text, out VascularNet r_vnet)
		{
            ILoadVascularNet vnet = new VascularNet();

            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CreateSpecificCulture("en-US");

            string[] readText = text.Split(new[] { '\r', '\n' });

            Regex regex = new Regex(@"^name:\s*(\w+)$", RegexOptions.IgnoreCase);
            int i = 0;
            while (!regex.IsMatch(readText[i]))
            {
                i++;
                if (i >= readText.Length)
                {
                   throw new TopReadingException("No correct name of vascular system in found\n");                   
                }
            }

            Match name_match = regex.Match(readText[i]);
            vnet.name = name_match.Groups[1].Value;

            

           
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


                        //regex = new Regex(@"^\s*(\d+)\s+X:(-*\d+.\d+)\s+Y:(-*\d+.\d+)\s+Z:(-*\d+.\d+)\s+R:(-*\d+.\d+)\s+C:(\d+.\d+)$", RegexOptions.IgnoreCase);
                        regex = new Regex(@"^\s*(\d+)\s+L:(-*\d+.\d+)\s+R:(-*\d+.\d+)\s+E:(-*\d+.\d+)\s+T:(-*\d+)$", RegexOptions.IgnoreCase);
                        Match node_match = regex.Match(readText[i]);

                        if (node_match.Groups.Count < 4)
                        {
                            //regex = new Regex(@"^\s*(\d+)\s+X:(-*\d+.\d+)\s+Y:(-*\d+.\d+)\s+Z:(-*\d+.\d+)\s+R:(-*\d+.\d+)$", RegexOptions.IgnoreCase);
                            regex = new Regex(@"^\s*(\d+)\s+L:(-*\d+.\d+)\s+R:(-*\d+.\d+)$", RegexOptions.IgnoreCase);
                            node_match = regex.Match(readText[i]);
                            if (node_match.Groups.Count < 4)
                            {
                                break;
                            }
                        }


                        int id = int.Parse(node_match.Groups[1].Value);
                        //convert to meter from millimeter
                        Vector1 position = new Vector1(double.Parse(node_match.Groups[2].Value)*1e-3);
                        double rad = double.Parse(node_match.Groups[3].Value);
                        double elasticity = double.Parse(node_match.Groups[4].Value) * 1e6;
                        int nodetype = int.Parse(node_match.Groups[5].Value);
                      
                        vnet.vascular_system.Add(new VascularNode(id, position, rad * 1e-3, elasticity, nodetype));
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
                            {
                                break;                                
                            }


                            int id = int.Parse(node_match[0].Value);
                            bonds_index.Add(new List<int>());
                            bonds_index[bonds_index.Count - 1].Add(id);

                            for (int n = 1; n < node_match.Count; n++)
                                bonds_index[bonds_index.Count - 1].Add(int.Parse(node_match[n].Value));

                            bond_string_count++;
                        }

                if (i >= readText.Length - 1)
                {                   
                    foreach (var str in bonds_index)
                    {
                        VascularNode nd = vnet.vascular_system.Find(x => x.id == str[0]);
                        for (int s = 1; s < str.Count; s++)
                            nd.addNeighbours(new VascularNode[] { vnet.vascular_system.Find(x => x.id == str[s]) });
                    }

                    break;
                }
            }
            r_vnet = (VascularNet)vnet;
            return false;
        }

        static public bool LoadTableFunctionFromString(string text, out TableFunction t_f)
        {
            t_f = null;

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo(localization);

            string[] readText = text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            if (readText.GetLength(0) == 0)
                return false;

            bool res = false;
            double[,] table_function_0 = new double[readText.Length, 2];
            double[,] table_function_1 = new double[readText.Length, 2];
            Regex regex = new Regex(@"^\s*(-?\d+(\.?\d*)?(e[-+]?\d+)?)$", RegexOptions.IgnoreCase);
            Match node_match = regex.Match(readText[0]);
            double scale = double.Parse(node_match.Groups[1].Value);
            for (int i = 1; i < table_function_1.GetLength(0); i++)
            {
                regex = new Regex(@"^\s*(-?\d+\.?\d*)[\s+\t]?(-?\d+\.?\d*)$", RegexOptions.IgnoreCase);
                node_match = regex.Match(readText[i]);
                double time = double.Parse(node_match.Groups[1].Value);
                double val = double.Parse(node_match.Groups[2].Value);
                table_function_1[i, 0] = time;
                table_function_1[i, 1] = val;
                table_function_0[i, 0] = time;
                table_function_0[i, 1] = val;
                res = true;
            };
            double period = table_function_0[table_function_0.GetLength(0)-1, 0];
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

        static public void WriteState(double current_time, IWriteVascularNet v_net, System.IO.StreamWriter out_dynamics_file)
        {
            string out_text = "";
            //Console.Write("Calc. time: ");
            //Console.WriteLine((current_time).ToString("F5"));
            
            out_text += "WT: ";
            out_text += current_time;
            out_text += "\n";

            for (int j = 0; j < v_net.vascular_system.Count; j++)
            {
                out_text += v_net.vascular_system[j].id.ToString();
                out_text += "\t";
                //out_text += (0.0f).ToString("F5");
                //out_text += "\t";
				out_text += (v_net.vascular_system[j].velocity * v_net.vascular_system[j].lumen_sq * 1e6).ToString("F8");
                out_text += "\t" + v_net.vascular_system[j].pressure.ToString("F8");
                out_text += "\t" + (Math.Sqrt(v_net.vascular_system[j].lumen_sq/Math.PI) * 1e3).ToString("F8");
                out_dynamics_file.WriteLine(out_text);
                out_text = "";
            }

            out_text += "\n\n";

            out_text = "";

            out_dynamics_file.Flush();         
        }

        static public void WriteSummaryState(double current_time, NodeSummary[] summary,  System.IO.StreamWriter out_summary_file)
        {
            string out_text = "";            
            
            out_text += "WT: ";
            out_text += current_time;
            out_text += "\n";

            for (int j = 0; j < summary.GetLength(0); j++)
            {
                out_text += summary[j].node.id.ToString();
                out_text += "\t";
                out_text += summary[j].proximal_dst.ToString();
                out_text += "\t";
                out_text += (summary[j].node.radius*1e3).ToString();
                out_text += "\t";                
                //out_text += (summary[j].get_av_flux()*1e6).ToString();
                //out_text += "\t";
                out_text += summary[j].Pd.ToString();
                out_text += "\t";
                out_text += summary[j].Ps.ToString();
                out_text += "\t";
                out_text += summary[j].get_av_prs().ToString();
                out_summary_file.WriteLine(out_text);
                out_text = "";

            }
            out_text += "\n\n";
            out_text = "";
            out_summary_file.Flush();         
        }

        static public TableFunction xScaleTableFunction(double timestep, double period, double new_period, TableFunction t_f)
        {
            int N = (int)Math.Ceiling(period / timestep);
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

        static public void WriteRCR(string filename, List<PressureOutletRCR> bc_list)
        {
            List<string> writeText = new List<string>();
            for (int i = 0; i < bc_list.Count; i++)
            {
                string out_text = "";
                out_text = out_text + bc_list[i].core_node.id.ToString() + " R1:" + (bc_list[i].R1 * 1e-9).ToString() + " R2:" + (bc_list[i].R2 * 1e-9).ToString() + " C:" + (bc_list[i].C * 1e12).ToString();
                writeText.Add(out_text);
            }
            File.WriteAllLines(filename, writeText);
        }

        static public bool LoadBC_RCR_paramsFromFile(string filename, out List<BC_params> BC_Params)
        {
            return LoadBC_RCR_paramsFromString(File.ReadAllText(filename), out BC_Params);
        }

        static public bool LoadBC_RCR_paramsFromString(string text, out List<BC_params> BC_Params)
        {
            BC_Params = new List<BC_params>();

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo(localization);
            string[] readText = text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

            bool res = false;
            Match node_match;
            for (int i = 0; i < readText.GetLength(0); i++)
            {
                //56          R1:  1.8104          R2:  7.2417         C:   31.29
                Regex regex = new Regex(@"^(\d+)\s+R1:(\d+\.\d+)\s+R2:(\d+\.\d+)\s+C:(\d+\.?\d*)$", RegexOptions.IgnoreCase);//^(\d+)[\s+\t]?R1:(\d+\.\d+)[\s+\t]?R2:(\d+\.\d+)[\s+\t]?C:(\d+\.\d+)$", RegexOptions.IgnoreCase);
                node_match = regex.Match(readText[i]);
                int id = int.Parse(node_match.Groups[1].Value);
                double R1 = double.Parse(node_match.Groups[2].Value);
                double R2 = double.Parse(node_match.Groups[3].Value);
                double C = double.Parse(node_match.Groups [4].Value);


                BC_params prms = new BC_params();
                BC_Params.Add(prms);
                
                prms.id = id;
                prms.R1 = R1 * 1e9;
                prms.R2 = R2 * 1e9;
                prms.C = C * 1E-12;

                BC_Params[i] = prms;

                res = true;
            };

            return res;
        }

        public static List<ClotDscr> readClotFile(string filename)
        {
            string text = File.ReadAllText(filename);
            string[] lines = text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            Regex regex = new Regex(@"^(\d+)\s+DG:(0.\d+)$", RegexOptions.IgnoreCase);
            List<ClotDscr> result = new List<ClotDscr>();
            for (int i = 0; i < lines.GetLength(0); i++)
            {
                Match clot_match = regex.Match(lines[i]);
                if (clot_match.Groups.Count != 3)
                    continue;

                int id = int.Parse(clot_match.Groups[1].Value);
                double DG = double.Parse(clot_match.Groups[2].Value);
                ClotDscr cd = new ClotDscr();
                cd.node_id = id;
                cd.degree = (float)DG;
                result.Add(cd);
            }

            return result;
        }

        public static void FlushStenosisDscr(string filename, List<ClotDscr> stenosis_list, double tot_seconds = 0)
        {
            List<string> out_strings = new List<string>();
            int curr_node_id = 0;
            for (int i = 0; i < stenosis_list.Count; i++ )
            {
                curr_node_id = stenosis_list[i].node_id;
                List<ClotDscr> sublist = stenosis_list.FindAll(x => x.node_id == curr_node_id);
                if (sublist.FindAll(x => !x.is_ready).Count != 0)
                    continue;

                string out_string = curr_node_id.ToString() + " ";
                out_string += sublist[0].reference_flux + " ";
                foreach (var s in sublist)
                {
                    out_string += s.degree + ":";
                    out_string += s.depressed_flux + " ";
                }
                out_strings.Add(out_string);
            }
            out_strings.Add("RunTime " + tot_seconds.ToString("0.0"));
            File.AppendAllLines(filename, out_strings);
            }

        public static void WriteStenosisDscr(string filename, List<ClotDscr> stenosis_list, double tot_seconds = 0)
        {
            List<string> out_strings = new List<string>();
            int curr_node_id = 0;;
            while (stenosis_list.Count > 0) 
            {
                curr_node_id = stenosis_list[0].node_id;
                List<ClotDscr> sublist = stenosis_list.FindAll(x => x.node_id == curr_node_id);               
                
                string out_string = curr_node_id.ToString()+ " ";
                out_string += sublist[0].reference_flux + " ";
                foreach (var s in sublist)
                {                    
                    out_string += s.degree + ":";
                    out_string += s.depressed_flux+ " ";
                    stenosis_list.Remove(s);
                }
                out_strings.Add(out_string);
                
            }
            out_strings.Add("RunTime " + tot_seconds.ToString("0.0"));
            File.WriteAllLines(filename, out_strings);            
        }
        public static void ReadClotFile(string filename, ref List<ClotNode> task)
        {
            string text = File.ReadAllText(filename);
            string[] clotnodes= text.Split(new[] {'\n' }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var st in clotnodes.Skip(1))
            {
                string[] clotinfo = st.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                ClotNode cn = new ClotNode();
                cn.Clot_id = int.Parse(clotinfo[0]);
                cn.Node_id = int.Parse(clotinfo[1]);
                cn.Permeability = float.Parse(clotinfo[2]);
                cn.Porosity= float.Parse(clotinfo[3]);
                task.Add(cn);
            }
        }

                public static void ReadModellingParameters(string filename)
        {
            string text = File.ReadAllText(filename);
            string[] lines = text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            Regex regex = new Regex(@"^(\w+)=\s*(.+)$", RegexOptions.IgnoreCase);
            // ugly but it works
            for (int line_index = 0; line_index < lines.GetLength(0); line_index++)
            {
                Match name_match = regex.Match(lines[line_index]);
                string name = name_match.Groups[1].Value;
                //double value = double.Parse(name_match.Groups[2].Value);

                if (name == "TIMESTEP")
                    Program.TIMESTEP = float.Parse(name_match.Groups[2].Value);

                if (name == "Beat_Duration")
                    Program.END_TIME = float.Parse(name_match.Groups[2].Value);

                if (name == "OUTPUT_PERIOD")
                    Program.OUTPUT_PERIOD = float.Parse(name_match.Groups[2].Value);
               
                if (name == "BLOOD_VISC")
                    BloodFlow.GlobalDefs.BLOOD_VISC = double.Parse(name_match.Groups[2].Value);

                //if (name == "YoungsModules")
                  //  BloodFlow.GlobalDefs.YOUNG_MODULUS = float.Parse(name_match.Groups[2].Value);

                if (name == "Density")
                    BloodFlow.GlobalDefs.BLOOD_DENSITY = double.Parse(name_match.Groups[2].Value);

                if (name == "DIASTOLIC_PRESSURE")
                    BloodFlow.GlobalDefs.DIASTOLIC_PRESSURE = double.Parse(name_match.Groups[2].Value);

                if (name == "OUT_PRESSURE")
                    BloodFlow.GlobalDefs.OUT_PRESSURE = double.Parse(name_match.Groups[2].Value);

                if (name == "FRICTION_C")
                    BloodFlow.GlobalDefs.FRICTION_C = double.Parse(name_match.Groups[2].Value);

                if (name == "ScriptLocation")
                    BloodFlow.GlobalDefs.ScriptLocation = name_match.Groups[2].Value;
            }
            Program.SimulationSteps = (int)Math.Ceiling(Program.END_TIME / Program.TIMESTEP);
            Program.TIMESTEP = Program.END_TIME / Program.SimulationSteps;
            //Program.OUTPUT_PERIOD = Program.OUTPUT_PERIOD * Program.TIMESTEP;
        }

        public static string readTaskFile(string filename, ref string top_filename, ref string par_filename, ref List<Tuple<int, string>> inlet_data, ref List<ClotDscr> task, ref string out_filename, ref int Coupled_Nodes_N)
        {
            string text = File.ReadAllText(filename);
            string[] lines = text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

			string[] path_tmp = filename.Split(new[] { Path.DirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries);

    //        string base_path = "";
    //        for (int i = 0; i < path_tmp.GetLength(0) - 1; i++)
				//base_path += path_tmp[i] + Path.DirectorySeparatorChar;

            int filenamelength = path_tmp.Last().Length;
            string base_path = filename.Substring(0, filename.Length - filenamelength);
            //modelling parameters
            string modelling_parameter_file = base_path + "Model_parameters.txt";
            ReadModellingParameters(modelling_parameter_file);

            Regex regex = new Regex(@"^<(\w+)>$", RegexOptions.IgnoreCase);

            int line_counter;
            for (line_counter = 0; line_counter < lines.GetLength(0); line_counter++)
            {
                Match name_match = regex.Match(lines[line_counter]);
                if (name_match.Groups.Count != 2)
                {                    
                    line_counter++;
                    continue;
                }
                out_filename = name_match.Groups[1].Value + ".out";
                line_counter++;
                break;
            }

            regex = new Regex(@"^Topology:\s*(\w+.*)$", RegexOptions.IgnoreCase);
            for (int i = line_counter; i < lines.GetLength(0); i++)
            {
                Match top_match = regex.Match(lines[i]);
                if (top_match.Groups.Count != 2)                
                    continue;
                
                top_filename = top_match.Groups[1].Value;
                break;
            }

            regex = new Regex(@"^InletFlux:\s+(.+)$", RegexOptions.IgnoreCase);
            for (int i = line_counter; i < lines.GetLength(0); i++)
            {
                Match inlet_match = regex.Match(lines[i]);
                if (inlet_match.Groups.Count != 2)                
                    continue;

                string[] inlets = inlet_match.Groups[1].Value.Split(',');
                foreach (var inlt in inlets)
                {   
                    string[] s = inlt.Trim().Split(' ');
                    Tuple<int, string> in_tlp = new Tuple<int, string>(int.Parse(s[0]), s[1].Trim());
                    inlet_data.Add(in_tlp);
                }

              //  inlet_filename = inlet_match.Groups[1].Value;
                break;
            }

            regex = new Regex(@"^OutletParams:\s*(\w+.*)$", RegexOptions.IgnoreCase);
            for (int i = line_counter; i < lines.GetLength(0); i++)
            {
                Match par_match = regex.Match(lines[i]);
                if (par_match.Groups.Count != 2)                
                    continue;
                
                par_filename = par_match.Groups[1].Value;
                break;
            }

            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo(localization);

            regex = new Regex(@"^Task:\s*(.+)$", RegexOptions.IgnoreCase);
            for (int i = line_counter; i < lines.GetLength(0); i++)
            {
                Match task_match = regex.Match(lines[i]);
                if (task_match.Groups.Count != 2)
                    continue;
                
                string[] stenosis =  task_match.Groups[1].Value.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                foreach(var st in stenosis)
                {
                   string [] id_dergee = st.Split(new[] { ':' }, StringSplitOptions.RemoveEmptyEntries);
                   ClotDscr dscr = new ClotDscr();
                   dscr.node_id =    int.Parse(id_dergee[0]);
                   dscr.degree  =  float.Parse(id_dergee[1]);
                   task.Add(dscr);                    
                }
            }

			regex = new Regex(@"^NumberCoupledPoints:\s*(.+)$", RegexOptions.IgnoreCase);
			for (int i = line_counter; i < lines.GetLength(0); i++)
            {
                Match par_match = regex.Match(lines[i]);
                if (par_match.Groups.Count != 2)
                    continue;

				//Console.WriteLine(par_match.Groups[1].Value);
				Coupled_Nodes_N=int.Parse(par_match.Groups[1].Value);
                break;
            }

            return base_path;
        }

        public static string localization;
    }
}
