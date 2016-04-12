package beast.app.seqgen;

import java.io.FileNotFoundException;
import java.io.PrintStream;

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
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.TreeTraitMap;//in beast-classic
import beast.util.Randomizer;
import beast.util.XMLProducer;
//import sphericalGeo.TreeTraitMap;

/**
 * Created by williamhsu on 18/03/16.
 */
@Description("Traits containing location data generated using " +
        "simple isotropic random walk are mapped onto a given tree.")
public class SimulatedTreeTraitMapBRW extends TreeTraitMap {

    //final public Input<TreeTraitMap> m_dataInput = new Input<>("data", "trait data which specifies... ", Validate.REQUIRED);
    //final public Input<Tree>m_treeInput = new Input<Tree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<String>m_traitNameInput = new Input<String>("traitName", "Name of trait to be map onto the tree");
    final public Input<Double>m_timeStepInput = new Input<Double>("timeStep", "time step between moves(default 0.01).", 0.01);
    final public Input<Double >m_spatialStepInput = new Input<Double>("spatialStep", "spatial step of move(default 0.01).", 0.01);
    public Input<Double> m_clockwiseProbInput = new Input<Double>("clockwiseProb", "probability of turning clock wise(default 0.5)", 0.4);
    public Input<Double> m_counterClockwiseProbInput = new Input<Double>("counterClockwiseProb", "probability of turning counter-clockwise(default 0.5)", 0.5);


    TreeInterface tree;
    RealParameter rP;
    String traitName;
    double [][] taxonLocations;
    double timeStep;
    double spatialStep;
    double direction;
    int numberOfLeaves;
    String value;
    double m_clockwiseProb;
    double m_counterClockwiseProb;

    double clockwiseCouterClockwise;

    public SimulatedTreeTraitMapBRW(){
        parameterInput.setRule(Validate.OPTIONAL);
    }

    public void initAndValidate(){
        tree = treeInput.get();
        rP = parameterInput.get();
        traitName = m_traitNameInput.get();
        timeStep = m_timeStepInput.get();
        spatialStep = m_spatialStepInput.get();
        m_clockwiseProb = m_clockwiseProbInput.get();
        m_counterClockwiseProb = m_counterClockwiseProbInput.get();

        if (!(tree instanceof Tree))
            throw new IllegalArgumentException("Tree input must be a true Tree, not just TreeInterface.");

        makeRootLocation((Tree)tree);//this needs to be a Tree object it is currently a tree interface
        value = "";
        numberOfLeaves = tree.getLeafNodeCount();
        for (int i = 0; i < numberOfLeaves; i++) {
            if (i == numberOfLeaves - 1) {
                value += "t" + Integer.toString(i) + "=" + Double.toString(taxonLocations[i][0]) + " " +
                        Double.toString(taxonLocations[i][1]) + "\n";
            } else {
                value += "t" + Integer.toString(i) + "=" + Double.toString(taxonLocations[i][0]) + " " +
                        Double.toString(taxonLocations[i][1]) + ", \n";
            }
        }
        setInputValue("value", value);
        super.initAndValidate();
    }

    private void makeRootLocation(Tree tree){
        Node root = tree.getRoot();
        System.out.println("Tree Height: " + root.getHeight());

        double [] rootLocation = {0.0, 0.0};

        taxonLocations = new double [tree.getNodeCount()][2];
        taxonLocations[root.getNr()] = rootLocation;

        //double angle = 0.0;
        BRW(root);
    }

    private void BRW(Node root){
        double rootHeight = root.getHeight();
        Node leftChild = root.getLeft();
        Node rightChild = root.getRight();
        double lTimeElapsed = rootHeight - leftChild.getHeight();
        double rTimeElapsed = rootHeight - rightChild.getHeight();

        //find number of steps and changes of directions on the left branch
        int numSteps;
        numSteps = (int)(lTimeElapsed/timeStep);
        //find location of left child
        double step;
        double lat;
        double lon;
        double cur_direction;
        double l_direction;
        //set left child location initially to that of parent
        taxonLocations[leftChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[leftChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            clockwiseCouterClockwise = Randomizer.nextDouble();
            cur_direction = Randomizer.nextDouble() * Math.PI;
            if (clockwiseCouterClockwise <= m_clockwiseProb)//direction is clockwise
                l_direction = cur_direction *-1;
            else
                l_direction = cur_direction;
            lat = Math.sin(l_direction)*spatialStep;
            lon = Math.cos(l_direction)*spatialStep;
            taxonLocations[leftChild.getNr()][0] = taxonLocations[leftChild.getNr()][0] + lat;
            taxonLocations[leftChild.getNr()][1] = taxonLocations[leftChild.getNr()][1] + lon;
        }
        //find number of steps and changes of directions on the right branch
        numSteps = (int)(rTimeElapsed/timeStep);
        double r_direction;
        //find location of right child
        //set right child location initially to that of parent
        taxonLocations[rightChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[rightChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            clockwiseCouterClockwise = Randomizer.nextDouble();
            cur_direction = Randomizer.nextDouble() *Math.PI;
            if (clockwiseCouterClockwise <= m_clockwiseProb)//direction is clockwise
                r_direction = cur_direction *-1;
            else
                r_direction = cur_direction;
            lat = Math.sin(r_direction)*spatialStep;
            lon = Math.cos(r_direction)*spatialStep;
            taxonLocations[rightChild.getNr()][0] = taxonLocations[rightChild.getNr()][0] + lat;
            taxonLocations[rightChild.getNr()][1] = taxonLocations[rightChild.getNr()][1] + lon;
        }
        if (leftChild.getChildCount() != 0){
            BRW(leftChild);
        }
        if (rightChild.getChildCount() != 0){
            BRW(rightChild);
        }
    }
}