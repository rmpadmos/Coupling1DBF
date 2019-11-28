using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using System.IO;

namespace BloodFlow
{
    public interface ILoadVascularNet
    {
        string name {get; set;}
        List<VascularNode> vascular_system {get; set;}
    }

    public interface IWriteVascularNet
    {
        string name { get; }
        List<VascularNode> vascular_system { get; }
    }

    public class VascularNet : ILoadVascularNet, IWriteVascularNet
    {
        static public int NumOfNeigbours(VascularNode node)
        {
            return node.neighbours.Count;
        }
        
        public double setHeartRate(int HR, double timestep)
        {
            if (HR > 120 || HR < 40)
                return 0;

            double new_period = 60.0 / HR;

            foreach (var bc in this.bounds)
                if (bc.GetType() == typeof(InletFlux))
                {
                    InletFlux inlet_bc = (InletFlux)bc; //TODO: replace by interface on inlet type of BCs
                    inlet_bc.flux_on_time = IO_Module.xScaleTableFunction(timestep, inlet_bc.base_period, new_period, inlet_bc.base_flux_on_time);
                }

            return new_period;
        }

        public VascularNet()
        {            
            vascular_system = new List<VascularNode>();

            knots = new List<Knot>();
            bounds = new List<BoundaryCondition>();
            threads = new List<Thread>();
        }

        public void defineNodeDirVectors(int[] initial_set, out getFloatValueDelegate getProximaDst, out setFloatValueDelegate setProximaDst)
        {
            float curr_value = 0;
            List<VascularNode> curr_front;
            curr_front = new List<VascularNode>();
            List<VascularNode> tmp_curr_front = new List<VascularNode>();

            foreach(var id in initial_set)
                curr_front.Add(vascular_system[0]);   
            
            getFloatValueDelegate get_del;
            setFloatValueDelegate set_del;
            VascularNode.newFloatValueLayer(out get_del, out set_del);
            foreach(var n in curr_front)
                set_del(n, curr_value);

            while (true)
            {
                foreach (var n in curr_front)
                    foreach (var ng in n.neighbours)
                    {
                        if (get_del(ng) > curr_value)
                        {
                            set_del(ng, curr_value);
                            tmp_curr_front.Add(ng);
                        }
                    }
                if (curr_front.Count == 0)
                    break;

                curr_front = new List<VascularNode>(tmp_curr_front);
                tmp_curr_front.Clear();
                curr_value+=1.0f;
            }

            foreach (var n in vascular_system)
            {
                if (n.neighbours.Count < 3)
                {
                    //Console.WriteLine(n.id);
                    //if (n.id == 1717)
                    //{
                    //    Console.WriteLine("1");
                    //}
                    //Console.WriteLine(get_del(n.neighbours.First()));
                    //Console.WriteLine(get_del(n));
                    if (n.neighbours.Last().position < n.neighbours.First().position)
                        n.neighbours.Reverse();
                    n.defDirVector();
                }
            }

            getProximaDst = get_del;
            setProximaDst = set_del;
        }

        public bool setCloth(int node_id, float degree)
        {            
            VascularNode nd = this.vascular_system.Find(x => x.id == node_id);

			try {
                Thread thr = this.node2thread[nd];

				foreach ( var child in nd.neighbours ) 
                {
					var th = this.node2thread[child];
				}
	            
				return(thr.setClot(nd, degree));
               // return (thr.setFFRClot(nd, degree));

			    } catch ( Exception ) {
				return false;
			}
        }

        public void removeClot(int node_id)
        {
            VascularNode nd = this.vascular_system.Find(x => x.id == node_id);
            Thread thr = this.node2thread[nd];
            thr.removeClot(nd);
        }      

        public List<VascularNode> getSubsystem(Predicate<VascularNode> p)
        {

            if (p == null)
                return vascular_system.GetRange(0, vascular_system.Count);
            else
            {
                List<VascularNode> subsystem = new List<VascularNode>();

                foreach (VascularNode n in vascular_system)
                    if (p(n))
                        subsystem.Add(n);
                return subsystem;
            }
        }
       
        public void defineNet(getFloatValueDelegate getProximalDistance, setFloatValueDelegate setProximalDistance, List<ClotNode> ClotNodes)
        {
            int count_id = vascular_system.Count;

            foreach (var n in vascular_system)
                n.defDirVector();


            // this we should get from import
            List<VascularNode> knot_nodes = getSubsystem(x => x.nodetype == 1);
            List<VascularNode> term_nodes = getSubsystem(x => x.nodetype == 2);
            for (int i = 0; i < ClotNodes.Count; i++)
            {
                ClotNodes[i].Node = vascular_system.Find(x => x.id == ClotNodes[i].Node_id);
            }
            foreach (var node in ClotNodes)
            {
                node.Node.nodetype = 3;
            }
            foreach (var node in vascular_system)
            {
                if (node.nodetype == 3)
                {
                    foreach (var nnode in node.neighbours)
                    {
                        if (nnode.nodetype <= 1)
                        {
                            node.nodetype = 4;
                            term_nodes.Add(node);
                        }
                    }
                }
            }

            // remove connections between clot boundary and the clot (terminal definition)
            foreach (var node in vascular_system)
            {
                if (node.nodetype >= 3)
                {
                    List<VascularNode> dellist = new List<VascularNode>();
                    foreach (var nnode in node.neighbours)
                    {   
                        if (nnode.nodetype == 4)
                        {
                            //term_nodes.Add(node);
                            dellist.Add(nnode);
                        }
                    }
                    foreach (var nnode in dellist)
                    {
                        node.excludeNeighbour(nnode);
                        nnode.excludeNeighbour(node);
                    }
                }
            }
            

            //List<int> list = new List<int>();
            //foreach (var node in ClotNodes)
            //{
            //    list.Add(node.Clot_id);
            //}
            //list = list.Distinct().ToList();

            ////List<ClotNode> clotboundary = new List<ClotNode>();
            //foreach (int clotid in list)
            //{
            //    List<ClotNode> listnodes = new List<ClotNode>();
            //    foreach (var node in ClotNodes)
            //    {
            //        if (node.Clot_id == clotid)
            //            listnodes.Add(node);
            //    }

            //    List<VascularNode> Node1 = getSubsystem(x => x.id == listnodes[0].Node_id);
            //    List<VascularNode> Node2 = getSubsystem(x => x.id == listnodes[listnodes.Count - 1].Node_id);

            //    //clotboundary.Add(listnodes[0]);
            //    //clotboundary.Add(listnodes[listnodes.Count - 1]);
            //}



            knots = new List<Knot>();
            bounds = new List<BoundaryCondition>();

            foreach (var tn in term_nodes)
            {
                bounds.Add(new BoundaryCondition(tn, 0));   // Just abstract BC, which do nothing             
            }

            foreach (var kn in knot_nodes)
            {knots.Add(new Knot(kn, 0)); // Just abstract knot, which do nothing and don't pass the solution through itself
            }

            //foreach (var kn in knot_nodes)
            //{
            //    List<VascularNode> b_nodes = new List<VascularNode>();
            //    List<VascularNode> n_nodes = new List<VascularNode>(kn.neighbours);

            //    // change this part to switch roles for temp nodes

            //    foreach (var n in n_nodes)
            //    {
            //        Vector1 pos = kn.position - (kn.position - n.position) * 0.001;
            //        double r = Math.Sqrt(n.lumen_sq_0 / Math.PI);
            //        VascularNode b_n = new VascularNode(count_id, pos, r);
            //        b_n.neighbours.Add(kn); 
            //        b_n.neighbours.Add(n);
            //        setProximalDistance(b_n, (getProximalDistance(kn) + getProximalDistance(n)) / 2);
            //        kn.neighbours.Remove(n);// do we need to remove neighbours?
            //        n.neighbours.Remove(kn);
            //        kn.addNeighbour(b_n);
            //        n.addNeighbour(b_n);
            //        b_nodes.Add(b_n);
            //        count_id++;
            //    }
            //    knots.Add(new Knot(kn, 0)); // Just abstract knot, which do nothing and don't pass the solution through itself
            //}

            List<string> writeText = new List<string>();

            threads     = new List<Thread>            ();
            node2thread = new Dictionary<VascularNode, Thread>();

            getBoolValueDelegate  isProcessed;
            setBoolValueDelegate setProcessed;
            VascularNode.newBoolValueLayer(out isProcessed, out setProcessed);

            VascularNode curr_node;
            //VascularNode next_node;

            foreach (var kn in knots)
                setProcessed(kn.core_node, true);        

            // Define threads starting from bounds
            foreach (var bc in bounds)
            {
                curr_node = bc.core_node;
                List<VascularNode> protothread = new List<VascularNode>();

                while(true)
                {
                    protothread.Add(curr_node);
                    setProcessed(curr_node, true);
                    try
                    {
                        curr_node = curr_node.neighbours.First(x => ( !isProcessed(x) ));
                    }
                    catch
                    {
                        break;
                    }
                }
                if (protothread.Count > 1)
                {
                    //if(getProximalDistance(protothread.First()) > getProximalDistance(protothread.Last()))
                    if (protothread[1].position > protothread.Last().position)
                    {
                        protothread.Reverse();
                    }
                    threads.Add(new Thread(protothread));
                    foreach (var n in protothread)                    
                        node2thread.Add(n, threads.Last());                    
                }
            }

            // Define threads starting from knots
            foreach (var kn in knots)
            {
                foreach (var start in kn.nodes)
                {
                    curr_node = start;                    
                    List<VascularNode> protothread = new List<VascularNode>();
                    while (true)
                    {
                        protothread.Add(curr_node);
                        setProcessed(curr_node, true);
                        try
                        {
                            curr_node = curr_node.neighbours.First(x => ( !isProcessed(x) ));
                        }
                        catch
                        {
                            break;
                        }
                    }
                    if (protothread.Count > 1)
                    {
                        //if(getProximalDistance(protothread.First()) > getProximalDistance(protothread.Last()))
                        if (protothread[1].position > protothread.Last().position)
                        {
                            protothread.Reverse();
                        }                
                        threads.Add(new Thread(protothread));
                        foreach (var n in protothread)
                            node2thread.Add(n, threads.Last());        
                    }
                }
            }            

            VascularNode.terminateBoolValueLayer(ref isProcessed);
        }

        public void specifyThreadType(int threadListid, Thread newTypeThread)
        {
            foreach (var n in this.threads[threadListid].nodes)            
                node2thread[n] = newTypeThread;

            this.threads[threadListid] = newTypeThread;                  
        }

      /*  public void specifyBoundaryType(BoundaryCondition oldTypeBC, BoundaryCondition newTypeBC)
        {
            foreach (var n in oldTypeBC.core_node)
                
                node2thread[n] = newTypeThread;
        }

        public void specifyKnotType(Thread oldTypeThread, Thread newTypeThread)
        {
            foreach (var n in oldTypeThread.nodes)
                node2thread[n] = newTypeThread;
        }*/

        public List<BoundaryCondition> getBounds()
        {
            return bounds;
        }

        public List<Thread> getThreads()
        {
            return threads;
        }
        
        public void fullReset()
        {
            foreach (var tr in threads)
            { tr.reset();}

            foreach (var kn in knots)
            { kn.reset();}

            foreach (var bc in bounds)
            { bc.reset();}
        }

        public string name {  get; set; }
        public List<VascularNode> vascular_system { get; set; }

        public List<Thread> threads;
        public List<Knot> knots;
        public List<BoundaryCondition> bounds;

        public Dictionary<VascularNode, Thread> node2thread;
    };
}
