package beast.brownianmotion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.List;

import beast.app.seqgen.SequenceSimulator;
import beast.core.Description;
import beast.core.Input;
import beast.core.BEASTObject;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.datatype.DataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import beast.util.XMLParser;
import beast.util.XMLProducer;

import java.awt.geom.Point2D;
import java.util.Dictionary;
import java.util.Enumeration;
import java.util.Hashtable;
/**
 * @author remco@cs.waikato.ac.nz
 */

/**
 * 
 * @author williamhsu
 * to run:  run configuration-> Name:  BrownianMotionSimulatorArray  Project:  MyBrownianMotion  Main class:  beast.app.beastapp.BeastMain
 * Arguments:  examples/testBrownianMotionSimulatorArray.xml
 */
@Description("Performs random sequence generation for a given site model. " +
		"Sequences for the leave nodes in the tree are returned as an alignment.")
public class BrownianMotionSimulatorArray extends beast.core.Runnable {
	public Input<RandomTree> random_treeInput = new Input<RandomTree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
	public Input<Double> varianceInput = new Input<Double>("variance", "variance of Gaussian distribution (default 1.0).", 1.0);
	public static Input<Integer>treesInput = new Input<Integer>("trees", "number of trees in this run (default 1).", 1);
	
	//put all node locations into a array
	static Point2D[] taxonLocations;
	
	private static void makeLocation(RandomTree tree, double variance){
		Node root = tree.getRoot();
		System.out.println(root.getHeight());
		Point2D rootLocation = new Point2D.Double(0,0);
		root.setLocation(rootLocation);
		taxonLocations = new Point2D[tree.getNodeCount()];
		taxonLocations[root.getNr()] = rootLocation;
		brownianMotion(root, variance);
	}
	
	
	private static void brownianMotion(Node root, double variance){
		Double rootHeight = root.getHeight();
		Node leftChild = root.getLeft();
		Node rightChild = root.getRight();
		Double lTimeElapsed = rootHeight - leftChild.getHeight();
		Double rTimeElapsed = rootHeight - rightChild.getHeight();


		double lx = Randomizer.nextGaussian()*Math.sqrt(lTimeElapsed*variance) + root.getLocation().getX();
		double ly = Randomizer.nextGaussian()*Math.sqrt(lTimeElapsed*variance) + root.getLocation().getY();			

		Point2D leftLocation = new Point2D.Double(lx,ly);
		leftChild.setLocation(leftLocation);
		
		taxonLocations[leftChild.getNr()]= leftLocation;
		
		//String taxonNrLeft = Integer.toString(leftChild.getNr());
		//System.out.println(leftChild.getNr() + "\t" +leftChild.getLocation());


		double rx = Randomizer.nextGaussian()*Math.sqrt(rTimeElapsed*variance) + root.getLocation().getX();
		double ry = Randomizer.nextGaussian()*Math.sqrt(rTimeElapsed*variance) + root.getLocation().getY();

		Point2D rightLocation = new Point2D.Double(rx,ry);
		rightChild.setLocation(rightLocation);
		
		taxonLocations[rightChild.getNr()]= rightLocation;
		
		//String taxonNrRight = Integer.toString(rightChild.getNr());
		//System.out.println(rightChild.getNr() + "\t" + rightChild.getLocation());

		if (leftChild.getChildCount() != 0){
			brownianMotion(leftChild, variance);
		}
		if (rightChild.getChildCount() != 0){
			brownianMotion(rightChild, variance);
		}

	}

	@Override
	public void initAndValidate() throws Exception {
		RandomTree randomTree = random_treeInput.get();
		double variance = varianceInput.get();
		int trees = treesInput.get();
	}



	@Override
	public void run() throws Exception {
		// TODO Auto-generated method stub
		//for(int i = 0; i< treesInput.get(); i++){
		makeLocation(random_treeInput.get(), varianceInput.get());
		System.out.println(random_treeInput.get());//print newick string
		System.out.println();
		for(Point2D x : taxonLocations){
			System.out.println(x);
		}
		
	
		//23/11/2015*************************************************
		//Enumeration<String>key = taxonLocationDict.keys();
		//while(key.hasMoreElements()){
		//	System.out.println(key.nextElement());
		//	}
		//Enumeration<Point2D>element = taxonLocationDict.elements();
		//while(element.hasMoreElements()){
		//	System.out.println(element.nextElement());
		//}
		//************************************************************
	}

	static TreeLikelihood getTreeLikelihood(BEASTObject plugin) throws Exception {
		for (BEASTObject plugin2 : plugin.listActivePlugins()) {
			if (plugin2 instanceof TreeLikelihood) {
				return (TreeLikelihood) plugin2;
			} else {
				TreeLikelihood likelihood = getTreeLikelihood(plugin2);
				if (likelihood != null) {
					return likelihood;
				}
			}
		}
		return null;
	}

	//
	public static void printUsageAndExit() {
		System.out.println("Usage: java " + SequenceSimulator.class.getName() + " <beast file> <nr of instantiations> [<output file>]");
		System.out.println("simulates from a treelikelihood specified in the beast file.");
		System.out.println("<beast file> is name of the path beast file containing the treelikelihood.");
		System.out.println("<nr of instantiations> is the number of instantiations to be replicated.");
		System.out.println("<output file> optional name of the file to write the sequence to. By default, the sequence is written to standard output.");
		System.exit(0);
	} // printUsageAndExit

	@SuppressWarnings("unchecked")
	public static void main(String[] args) {
		try {
			// parse arguments
			if (args.length < 2) {
				printUsageAndExit();
			}
			String sFile = args[0];
			//int nReplications = Integer.parseInt(args[1]);
			PrintStream out = System.out;
			if (args.length == 3) {
				File file = new File(args[2]);
				out = new PrintStream(file);
			}

			// grab the file
			String sXML = "";
			BufferedReader fin = new BufferedReader(new FileReader(sFile));
			while (fin.ready()) {
				sXML += fin.readLine();
			}
			fin.close();

			// parse the xml
			XMLParser parser = new XMLParser();
			BEASTObject plugin = parser.parseFragment(sXML, true);

			//find relevant objects from the model
			TreeLikelihood treeLikelihood = getTreeLikelihood(plugin);
			if (treeLikelihood == null) {
				throw new Exception("No treelikelihood found in file. Giving up now.");
			}
			//Alignment data = ((Input<Alignment>) treeLikelihood.getInput("data")).get();
			Tree tree = ((Input<Tree>) treeLikelihood.getInput("tree")).get();
			//SiteModel pSiteModel = ((Input<SiteModel>) treeLikelihood.getInput("siteModel")).get();
			//BranchRateModel pBranchRateModel = ((Input<BranchRateModel>) treeLikelihood.getInput("branchRateModel")).get();


			// feed to sequence simulator and generate leaves
				BrownianMotionSimulatorArray bm_Simulator = new BrownianMotionSimulatorArray();
				bm_Simulator.init(tree);
			//XMLProducer producer = new XMLProducer();
			//Alignment alignment = treeSimulator.simulate();
			//sXML = producer.toRawXML(alignment);
			//out.println("<beast version='2.0'>");
			//out.println(sXML);
			//out.println("</beast>");
		} catch (Exception e) {
			e.printStackTrace();
		}
	} // main

} // class SequenceSimulator

