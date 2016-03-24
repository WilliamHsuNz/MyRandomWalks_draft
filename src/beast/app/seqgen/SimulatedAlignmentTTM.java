package beast.app.seqgen;

/**
 * Created by williamhsu on 18/03/16.
 */

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.XMLProducer;
import sphericalGeo.TreeTraitMap;


/**
 * @author remco@cs.waikato.ac.nz
 */
@Description("An alignment containing sequences randomly generated using a"
        + "given site model down a given tree.")
public class SimulatedAlignmentTTM extends Alignment {
    final public Input<Alignment> m_data = new Input<>("data", "alignment data which specifies datatype and taxa of the beast.tree", Validate.REQUIRED);
    final public Input<Tree> m_treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<SiteModel.Base> m_pSiteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> m_pBranchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Integer> m_sequenceLengthInput = new Input<>("sequencelength", "nr of samples to generate (default 1000).", 1000);
    final public Input<String> m_outputFileNameInput = new Input<>(
            "outputFileName",
            "If provided, simulated alignment is additionally written to this file.");
    final public Input<Double>m_timeStepInput = new Input<Double>("timeStep", "time step of random walk (default 0.01).", 0.01);
    final public Input<Double>m_spatialStepInput = new Input<Double>("spatialStep", "spatial step of random walk (default 0.01).", 0.01);
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
    double [][]taxonLocations;
    double timeStep;
    double spatialStep;
    double direction;

    /**
     * an array used to transfer transition probabilities
     */
    protected double[][] m_probabilities;

    public SimulatedAlignmentTTM() {

        // Override the sequence input requirement.
        sequenceInput.setRule(Validate.OPTIONAL);
    }

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
        timeStep = m_timeStepInput.get();
        spatialStep = m_spatialStepInput.get();

        sequenceInput.get().clear();

        simulate();

        // Write simulated alignment to disk if requested:
        if (m_outputFileName != null) {
            PrintStream pstream;
            try {
                pstream = new PrintStream(m_outputFileName);
                pstream.println(new XMLProducer().toRawXML(this));
                pstream.close();
            } catch (FileNotFoundException e) {
                throw new IllegalArgumentException(e.getMessage());
            }
        }

        super.initAndValidate();
    }

    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    Sequence intArray2Sequence(int[] seq, Node node) {
        DataType dataType = m_data.get().getDataType();
        String seqString = dataType.state2string(seq);
//    	StringBuilder seq = new StringBuilder();
//    	String map = m_data.get().getMap();
//    	if (map != null) {
//    		for (int i  = 0; i < m_sequenceLength; i++) {
//    			seq.append(map.charAt(seq[i]));
//    		}
//    	} else {
//    		for (int i  = 0; i < m_sequenceLength-1; i++) {
//    			seq.append(seq[i] + ",");
//    		}
//			seq.append(seq[m_sequenceLength-1] + "");
//    	}
        String taxon = m_data.get().getTaxaNames().get(node.getNr());
        return new Sequence(taxon, seqString);
    } // intArray2Sequence

    /**
     * perform the actual sequence generation
     *
     * @return alignment containing randomly generated sequences for the nodes in the
     *         leaves of the tree
     */
    public void simulate() {
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

        //alignment.setDataType(m_siteModel.getFrequencyModel().getDataType());

        traverse(root, seq, category);
        geoSpread(m_tree);
        for(double [] x : taxonLocations)
            System.out.println(Arrays.toString(x));

        Double [] ttm_parameters = new Double [taxonLocations.length * taxonLocations[0].length];//taxonLocations.length
        for(int i = 0; i < taxonLocations.length; i++){
            //lat at the start 1/2 of the array
            ttm_parameters[i] = taxonLocations[i][0];
            //lon at the end 1/2 of the array
            ttm_parameters[i + taxonLocations.length] = taxonLocations[i][1];
        }

        RealParameter rp = new RealParameter(ttm_parameters);
        TreeTraitMap ttm =  new TreeTraitMap();
        ttm.init(m_tree, rp, "location");
        ttm.setInputValue("randomizelower", "-90, -180");
        ttm.setInputValue("randomizeupper", "90, 180");
        System.out.println("Leaf count: " + m_tree.getLeafNodeCount());
        String value = "";
        for(int i = 0; i < m_tree.getLeafNodeCount(); i++){
            value += "t"+ Integer.toString(i)+"_"+ Double.toString(taxonLocations[i][0])+"_"+
            Double.toString(taxonLocations[i][1])+ "=" + Double.toString(taxonLocations[i][0])+ " "+
            Double.toString(taxonLocations[i][1])+", ";
        }
        //System.out.println(value);
        ttm.setInputValue("value", value);

        System.out.println("value: " + ttm.value.get());

    } // simulate

    private void geoSpread(Tree tree){
        Node root = tree.getRoot();
        System.out.println("Tree Height: " + root.getHeight());
        //System.out.println("Random Walk Type " + m_rwType);

        double [] rootLocation = {0.0, 0.0};

        taxonLocations = new double [tree.getNodeCount()][2];
        taxonLocations[root.getNr()] = rootLocation;

        //double angle = 0.0;
        RW(root);
    }

    private void RW(Node root){
        double rootHeight = root.getHeight();
        Node leftChild = root.getLeft();
        Node rightChild = root.getRight();
        double lTimeElapsed = rootHeight - leftChild.getHeight();
        double rTimeElapsed = rootHeight - rightChild.getHeight();

        //find number of steps and changes of directions on the left branch
        int numSteps = 0;
        numSteps = (int)(lTimeElapsed/timeStep);
        //find location of left child
        double step;
        double lat;
        double lon;
        //set left child location initially to that of parent
        taxonLocations[leftChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[leftChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            direction = Randomizer.nextDouble() *Math.PI;
            direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counter-clockwise
            //step = l_steps[i];
            lat = Math.sin(direction)*spatialStep;
            lon = Math.cos(direction)*spatialStep;
            taxonLocations[leftChild.getNr()][0] = taxonLocations[leftChild.getNr()][0] + lat;
            taxonLocations[leftChild.getNr()][1] = taxonLocations[leftChild.getNr()][1] + lon;
        }
        //find number of steps and changes of directions on the right branch
        numSteps = (int)(rTimeElapsed/timeStep);
        //find location of right child
        //set right child location initially to that of parent
        taxonLocations[rightChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[rightChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            direction = Randomizer.nextDouble() *Math.PI;
            direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counterclockwise
            lat = Math.sin(direction)*spatialStep;
            lon = Math.cos(direction)*spatialStep;
            taxonLocations[rightChild.getNr()][0] = taxonLocations[rightChild.getNr()][0] + lat;
            taxonLocations[rightChild.getNr()][1] = taxonLocations[rightChild.getNr()][1] + lon;
        }
        if (leftChild.getChildCount() != 0){
            RW(leftChild);
        }
        if (rightChild.getChildCount() != 0){
            RW(rightChild);
        }
    }


    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     */
    void traverse(Node node, int[] parentSequence, int[] category) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
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
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category);
            }
        }
    } // traverse

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

//        // TODO Hack until SiteRateModel issue is resolved
//        if (m_siteModel.getSubstitutionModel() instanceof SubstitutionEpochModel) {
//            ((SubstitutionEpochModel)m_siteModel.getSubstitutionModel()).getTransitionProbabilities(tree.getNodeHeight(node),
//                    tree.getNodeHeight(parent),branchLength, probs);
//            return;
//        }
        //m_siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probs);
        m_siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);

    } // getTransitionProbabilities


} // class SequenceAlignment

