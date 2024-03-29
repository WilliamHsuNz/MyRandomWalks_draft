package beast.app.seqgen;

//import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.List;
import java.util.Arrays;

import beast.core.BEASTInterface;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

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


	
/**
 * @author William Hsu
 * Extension on seqgen by remco@cs.waikato.ac.nz to simulate geo spread based on (R)RW brownian motion
 * 15/2/2015
 * The programme simulates location data by making a (R)RW as outlined by Lemey2010
 * Random samples from a gamma/log-normal distribution are used as scalers for the covariance matrix 
 * in simulating RRW
 * Locations are then drawn from the Multivariate Normal distribution
 * The class produces a Nexus file of simulated sequence and location. 
 * 
 * To run:  Run configurations -> Name: seqgenRRWLogNormal  Project: MyBrownianMotion  
 * Main class: beast.app.beastapp.BeastMain  
 * Arguments:  examples/testSeqGenRRWLogNormal.xml
 * 
 */
@Description("Performs random sequence generation for a given site model. " +
        "Sequences for the leave nodes in the tree are returned as an alignment.")
public class SeqgenRWGammaLogNormal extends beast.core.Runnable {
    public Input<Alignment> m_data = new Input<Alignment>("data", "alignment data which specifies datatype and taxa of the beast.tree", Validate.REQUIRED);
    public Input<Tree> m_treeInput = new Input<Tree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    public Input<SiteModel.Base> m_pSiteModelInput = new Input<SiteModel.Base>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<BranchRateModel.Base>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    public Input<Integer> m_sequenceLengthInput = new Input<Integer>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    public Input<String> m_outputFileNameInput = new Input<String>("outputFileName","If provided, simulated alignment is written to this file rather " + "than to standard out.");
	
    public Input<String> m_rwTypeInput = new Input<String>("rwType", "Random Walk tyep");
    public Input<Double> v1_varianceInput = new Input<Double>("v1_variance", "variance of Normal distribution in 1st dimension (default 1.0).", 1.0);
    public Input<Double> v2_varianceInput = new Input<Double>("v2_variance", "variance of Normal distribution in 2nd dimension (default 1.0).", 1.0);
    public Input<Double> correlationCoefficientInput = new Input<Double>("correlationCoefficient", "correlation for specifying the covariance matrix (default 0.0.", 0.0);
    public Input<Double> scaleInput = new Input<Double>("scale", "scale of log-normal distribution (default 0.0) ", 0.0);
    public Input<Double> shapeInput = new Input<Double>("shape", "shape of log-normal distribution (default 1.0", 1.0);
    
    /**
     * nr of samples to generate *
     */
    protected int m_sequenceLength;
    /**
     * tree used for generating samples *
     */
    protected Tree m_tree;
    /**
     * site model used for generating samples *
     */
    protected SiteModel.Base m_siteModel;
    /**
     * branch rate model used for generating samples *
     */
    protected BranchRateModel m_branchRateModel;
    /**
     * nr of categories in site model *
     */
    int m_categoryCount;
    /**
     * nr of states in site model *
     */
    int m_stateCount;
    
    /**
     * name of output file *
     */
    String m_outputFileName;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[][] m_probabilities;
    
    /***
     * string defining type of RW model
     */
    static String m_rwType;
    
    /**
     *  variance of Brownian Motion
     */
    static double v1_variance;
    static double v2_variance;
    static double correlationCoefficient;

    /**
     * parameters for log-normal distribution
     */
    static double scale;
    static double shape;
    
    public static double [][] taxonLocations;
    static double [][] sigma; 
   

    @Override
    public void initAndValidate() {
        m_tree = m_treeInput.get();
        m_siteModel = m_pSiteModelInput.get();
        m_branchRateModel = m_pBranchRateModelInput.get();
        m_sequenceLength = m_sequenceLengthInput.get();
        m_stateCount = m_data.get().getMaxStateCount();
        m_categoryCount = m_siteModel.getCategoryCount();
        m_probabilities = new double[m_categoryCount][m_stateCount * m_stateCount];
        m_outputFileName = m_outputFileNameInput.get();
        m_rwType = m_rwTypeInput.get();
		v1_variance = v1_varianceInput.get();
		v2_variance = v2_varianceInput.get();
		correlationCoefficient = correlationCoefficientInput.get();
		scale = scaleInput.get();
		shape = shapeInput.get();
    }
    
    
    @Override	
    public void run() throws Exception {
    	//System.out.println(m_treeInput.get());//print newick String
    	makeLocation(m_tree);
    	System.out.println(m_tree);
    	
		for(double [] x : taxonLocations)
			System.out.println(Arrays.toString(x));
		
        Alignment alignment = simulate();
        
        //Write output to stdout or file
        PrintStream pstream;
        if (m_outputFileName == null)
            pstream = System.out;
        else
        	pstream = new PrintStream(m_outputFileName);        	
        pstream.println(new XMLProducer().toRawXML(alignment));
        
        //write Nexus file using m_outputFileName
        File file = new File("/Users/williamhsu/Documents/workspace/MyBrownianMotion/seqgenRW"+ m_rwType+ ".nex");
        file.getParentFile().mkdirs();
        PrintWriter writer = new PrintWriter(file);
    	writer.println("#NEXUS");
    	writer.println("BEGIN DATA;");
    	writer.println("       DIMENSIONS  NTAX =" + m_tree.getLeafNodeCount() + " NCHAR=" + m_sequenceLength + ";");
    	writer.println("       FORMAT DATATYPE = DNA  GAP = - MISSING = ?;");
    	writer.println("       MATRIX");
    	
    	BufferedReader br = new BufferedReader(new FileReader(m_outputFileName));
    	try{
    		String line = br.readLine();
    		while(line != null){
    			if (line.matches("<data>")||line.matches("</data>")){
    				line = br.readLine();
    			}else{
    				String delims = "[']";
    				String[] nexusValues = line.split(delims);
    				if (nexusValues.length > 4){
    					writer.println(nexusValues[1]);
    					writer.println(nexusValues[3]);
    				}	
    				line = br.readLine();
    			} 			
    		}	 
    	}finally{
    		br.close();
    	}
    	writer.println(";");
    	writer.println("End");
    	writer.close();
    }
    
    /**
     *Brownian Motion
     */	
    
	private static void makeLocation(Tree tree){
		Node root = tree.getRoot();
		System.out.println("Tree Height: " + root.getHeight());
		System.out.println("Random Walk Type " + m_rwType);
		System.out.println("v1_Variance: " + v1_variance);
		System.out.println("v2_Variance: " + v2_variance);
		System.out.println("CorrelationCoefficient: " + correlationCoefficient);
		System.out.println("Scale: " + scale);
		System.out.println("Shape: " + shape);
		
		double [] rootLocation = {0.0, 0.0};
		double [][] sigma = {{v1_variance, correlationCoefficient},{correlationCoefficient, v2_variance}};
		
		taxonLocations = new double [tree.getNodeCount()][];
		taxonLocations[root.getNr()] = rootLocation;
		
		brownianMotion(root, sigma);
	}
	
	private static void brownianMotion(Node root, double [][] sigma){
		double rootHeight = root.getHeight();
		Node leftChild = root.getLeft();
		Node rightChild = root.getRight();
		double lTimeElapsed = rootHeight - leftChild.getHeight();
		double rTimeElapsed = rootHeight - rightChild.getHeight();
		double l_d_sample;
		double r_d_sample;
		
		if(m_rwType == "Gamma"){
			//set up a gamma distribution for RRW
			GammaDistribution gd = new GammaDistribution(shape, scale);
			//draw a random sample from the gamma distribution
			l_d_sample = gd.sample();
		}if(m_rwType == "LogNormal"){
			//set up a log normal distribution for RRW
			LogNormalDistribution lnd = new LogNormalDistribution(scale, shape);
			//draw a random sample from the log normal distribution
			l_d_sample = lnd.sample();
		}else{
			//set up SRW
			l_d_sample = 1.0;
		}
		//scale sigma in time and l_d_sample
		double [][]l_sigma = new double [sigma.length][sigma[0].length];
		double l_RW_scaler = lTimeElapsed * l_d_sample;  
		for (int i = 0; i < sigma.length; i++){
			for (int j = 0; j < sigma[i].length; j++){
				l_sigma[i][j] = l_RW_scaler * sigma[i][j];
			}
		}		
		//sampling from multivariate normal distribution
		MultivariateNormalDistribution l_mnd = new MultivariateNormalDistribution(taxonLocations[root.getNr()], l_sigma);
		taxonLocations[leftChild.getNr()] = l_mnd.sample();	
		//change nr to include location
		//System.out.println(leftChild.getNr());
		
		if(m_rwType == "Gamma"){
			//set up a gamma distribution for RRW
			GammaDistribution gd = new GammaDistribution(shape, scale);
			//draw a random sample from the gamma distribution
			r_d_sample = gd.sample();
		}if(m_rwType == "LogNormal"){
			//set up a log normal distribution for RRW
			LogNormalDistribution lnd = new LogNormalDistribution(scale, shape);
			//draw a random sample from the log normal distribution
			r_d_sample = lnd.sample();
		}else{
			//set up SRW
			r_d_sample = 1.0;
		}
		//scale sigma in time and r_lnd_sample();
		double [][]r_sigma = new double [sigma.length][sigma[0].length];
		double r_RW_scaler = rTimeElapsed * r_d_sample;
		for (int i = 0; i < sigma.length; i++){
			for (int j = 0; j < sigma[i].length; j++){
				r_sigma[i][j] = r_RW_scaler * sigma[i][j];
			}
		}		
		//sampling from multivariate normal distribution
		MultivariateNormalDistribution r_mnd = new MultivariateNormalDistribution(taxonLocations[root.getNr()], r_sigma);
		taxonLocations[rightChild.getNr()] = r_mnd.sample();

		if (leftChild.getChildCount() != 0){
			brownianMotion(leftChild, sigma);
		}
		if (rightChild.getChildCount() != 0){
			brownianMotion(rightChild, sigma);
		}

	}


    /**
     * Convert integer representation of sequence into a Sequence
     *`
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     * @throws Exception
     */
    Sequence intArray2Sequence(int[] seq, Node node) throws Exception {
        DataType dataType = m_data.get().getDataType();
        String sSeq = dataType.state2string(seq);
//    	StringBuilder sSeq = new StringBuilder();
//    	String sMap = m_data.get().getMap();
//    	if (sMap != null) {
//    		for (int i  = 0; i < m_sequenceLength; i++) {
//    			sSeq.append(sMap.charAt(seq[i]));
//    		}
//    	} else {
//    		for (int i  = 0; i < m_sequenceLength-1; i++) {
//    			sSeq.append(seq[i] + ",");
//    		}
//			sSeq.append(seq[m_sequenceLength-1] + "");
//    	}
        List<Sequence> taxa = m_data.get().sequenceInput.get();
        
        String sTaxon = taxa.get(node.getNr()).taxonInput.get();
        
        //attaching location x and y at the end of the sTaxon
        double []location = taxonLocations[Integer.parseInt(sTaxon.replace("t", ""))];
        double x = location[0];
        double y = location[1];
        sTaxon = sTaxon + "_" + x + "_" + y;
 
        return new Sequence(sTaxon, sSeq.toString());
    } // intArray2Sequence

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     * @throws Exception
     */
    public Alignment simulate() throws Exception {
        Node root = m_tree.getRoot();

        double[] categoryProbs = m_siteModel.getCategoryProportions(root);
        int[] category = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        double[] frequencies = m_siteModel.getSubstitutionModel().getFrequencies();
        int[] seq = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            seq[i] = Randomizer.randomChoicePDF(frequencies);
        }
        
        Alignment alignment = new Alignment();
        //alignment.setDataType(m_siteModel.getFrequencyModel().getDataType());
        //alignment.userDataTypeInput.setValue(m_data.get().getDataType(), alignment);
        //alignment.setID("SequenceSimulatorBrownianMotion");
        traverse(root, seq, category, alignment);
      
        return alignment;
    } 

    
    
    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     * @throws Exception
     */
    void traverse(Node node, int[] parentSequence, int[] category, Alignment alignment) throws Exception {
        for (int iChild = 0; iChild < 2; iChild++) {
            Node child = (iChild == 0 ? node.getLeft() : node.getRight());
            for (int i = 0; i < m_categoryCount; i++) {
                getTransitionProbabilities(m_tree, child, i, m_probabilities[i]);
            }

            int[] seq = new int[m_sequenceLength];
            double[] cProb = new double[m_stateCount];
            for (int i = 0; i < m_sequenceLength; i++) {
                System.arraycopy(m_probabilities[category[i]], parentSequence[i] * m_stateCount, cProb, 0, m_stateCount);
                seq[i] = Randomizer.randomChoicePDF(cProb);
            }

            if (child.isLeaf()) {
                alignment.sequenceInput.setValue(intArray2Sequence(seq, child), alignment);
            } else {
                traverse(child, seq, category, alignment);
            }
        }
    }

    /**
     * get transition probability matrix for particular rate category *
     */
    void getTransitionProbabilities(Tree tree, Node node, int rateCategory, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(node));
        branchRate *= m_siteModel.getRateForCategory(rateCategory, node);

        // Get the operational time of the branch

        //final double branchTime = branchRate * (parent.getHeight() - node.getHeight());

        //if (branchTime < 0.0) {
        //    throw new RuntimeException("Negative branch length: " + branchTime);
        //}

        //double branchLength = m_siteModel.getRateForCategory(rateCategory) * branchTime;

//        //TODO Hack until SiteRateModel issue is resolved
//        if (m_siteModel.getSubstitutionModel() instanceof SubstitutionEpochModel) {
//            ((SubstitutionEpochModel)m_siteModel.getSubstitutionModel()).getTransitionProbabilities(tree.getNodeHeight(node),
//                    tree.getNodeHeight(parent),branchLength, probs);
//            return;
//        }
        //m_siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probs);
        m_siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);

    } // getTransitionProbabilities
  
    /**
     * find a treelikelihood object among the plug-ins by recursively inspecting plug-ins *
     */
//    static TreeLikelihood getTreeLikelihood(BEASTObject plugin) throws Exception {
//        for (BEASTObject plugin2 : plugin.listActivePlugins()) {
//            if (plugin2 instanceof TreeLikelihood) {
//                return (TreeLikelihood) plugin2;
//            } else {
//                TreeLikelihood likelihood = getTreeLikelihood(plugin2);
//                if (likelihood != null) {
//                    return likelihood;
//                }
//            }
//        }
//        return null;
//    }
    static TreeLikelihood getTreeLikelihood(BEASTInterface plugin) throws Exception {
        for (BEASTInterface plugin2 : plugin.listActiveBEASTObjects()) {
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

    /**
     * helper method *
     */
    public static void printUsageAndExit() {
        System.out.println("Usage: java " + SeqgenRWGammaLogNormal.class.getName() + " <beast file> <nr of instantiations> [<output file>]");
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
            int nReplications = Integer.parseInt(args[1]);
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
            BEASTInterface plugin = parser.parseFragment(sXML, true);

            // find relevant objects from the model
            TreeLikelihood treeLikelihood = getTreeLikelihood(plugin);
            if (treeLikelihood == null) {
                throw new Exception("No treelikelihood found in file. Giving up now.");
            }
            Alignment data = ((Input<Alignment>) treeLikelihood.getInput("data")).get();
            Tree tree = ((Input<Tree>) treeLikelihood.getInput("tree")).get();
            SiteModel pSiteModel = ((Input<SiteModel>) treeLikelihood.getInput("siteModel")).get();
            BranchRateModel pBranchRateModel = ((Input<BranchRateModel>) treeLikelihood.getInput("branchRateModel")).get();


            // feed to sequence simulator and generate leaves
            SeqgenRWGammaLogNormal treeSimulator = new SeqgenRWGammaLogNormal();
            treeSimulator.init(data, tree, pSiteModel, pBranchRateModel, nReplications);
            XMLProducer producer = new XMLProducer();
            Alignment alignment = treeSimulator.simulate();
            sXML = producer.toRawXML(alignment);
            out.println("<beast version='2.0'>");
            out.println(sXML);
            out.println("</beast>");
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    } // main

} // class SequenceSimulator

