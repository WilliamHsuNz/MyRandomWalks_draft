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
 * This Programme extends on The TreeTraitMap Class from the beast-classic package.
 * The programme simulates Simple Isotropic Random Walk with boundaries down a
 * phylogenetic tree and maps geographical location traits to the nodes of the tree.
 * The programme then performs MCMC sampliing method to infer the location of the root
 * from the locations of the leaves.
 *
 * The programme can is run through the main method at beast.app.beastapp.BeastMain
 * The programme takes examples/testSimulatedTreeTraitMapWithBoundaries.xml as input arguement
 */
@Description("Traits containing location data generated using " +
        "simple isotropic random walk are mapped onto a given tree.")
public class SimulatedTreeTraitMapWithBoundaries extends TreeTraitMap {

    final public Input<String>m_traitNameInput = new Input<String>("traitName", "Name of trait to be map onto the tree");
    final public Input<Double>m_timeStepInput = new Input<Double>("timeStep", "time step between moves(default 0.01).", 0.01);
    final public Input<Double >m_spatialStepInput = new Input<Double>("spatialStep", "spatial step of move(default 0.01).", 0.01);
    public Input<Double> m_maxLatBoundaryInput = new Input<Double>("maxLatBoundary", "Maximum Latitude boundary (default 20.0).",1.0);
    public Input<Double> m_minLatBoundaryInput = new Input<Double>("minLatBoundary", "Maximum Latitude boundary (default -20.0).", -1.0);
    public Input<Double> m_maxLonBoundaryInput = new Input<Double>("maxLonBoundary", "Maximum Longitude boundary(default 20.0).", 1.0);
    public Input<Double> m_minLonBoundaryInput = new Input<Double>("minLonBoundary", "Minimum Longitude boundary(default -20.0).", -1.0);

    TreeInterface tree;
    RealParameter rP;
    String traitName;
    double [][] taxonLocations;
    double timeStep;
    double spatialStep;
    double direction;
    int numberOfLeaves;
    String value;
    double maxLatBoundary;
    double minLatBoundary;
    double maxLonBoundary;
    double minLonBoundary;

    public SimulatedTreeTraitMapWithBoundaries(){
        parameterInput.setRule(Validate.OPTIONAL);
    }

    public void initAndValidate(){
        tree = treeInput.get();
        rP = parameterInput.get();
        traitName = m_traitNameInput.get();
        timeStep = m_timeStepInput.get();
        spatialStep = m_spatialStepInput.get();
        maxLatBoundary = m_maxLatBoundaryInput.get();
        minLatBoundary = m_minLatBoundaryInput.get();
        maxLonBoundary = m_maxLonBoundaryInput.get();
        minLonBoundary = m_minLonBoundaryInput.get();

        if (!(tree instanceof Tree))
            throw new IllegalArgumentException("Tree input must be a true Tree, not just TreeInterface.");

        Node root = tree.getRoot();
        double [] rootLocation = {0.0, 0.0};
        taxonLocations = new double[tree.getNodeCount()][2];
        taxonLocations[root.getNr()] = rootLocation;
        SRWWithBoundaries(root);

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

    private void SRWWithBoundaries(Node root){
        double rootHeight = root.getHeight();
        Node leftChild = root.getLeft();
        Node rightChild = root.getRight();
        double lTimeElapsed = rootHeight - leftChild.getHeight();
        double rTimeElapsed = rootHeight - rightChild.getHeight();

        //find number of steps and changes of directions on the left branch
        int numSteps = 0;
        numSteps = (int)(lTimeElapsed/timeStep);
        //find location of left child
        //double step;
        double lat;
        double lon;
        double proposedLat;
        double proposedLon;
        //set left child location initially to that of parent
        taxonLocations[leftChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[leftChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            direction = Randomizer.nextDouble() * Math.PI;
            direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counter-clockwise
            lat = Math.sin(direction)*spatialStep;
            lon = Math.cos(direction)*spatialStep;
            proposedLat = taxonLocations[leftChild.getNr()][0] + lat;
            proposedLon = taxonLocations[leftChild.getNr()][1] + lon;
            //boundary handling
            while (proposedLat > maxLatBoundary || proposedLat < minLatBoundary || proposedLon > maxLonBoundary || proposedLon < minLonBoundary){
                direction = Randomizer.nextDouble() *Math.PI;
                direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counter-clockwise
                lat = Math.sin(direction)*spatialStep;
                lon = Math.cos(direction)*spatialStep;
                proposedLat = taxonLocations[leftChild.getNr()][0] + lat;
                proposedLon = taxonLocations[leftChild.getNr()][1] + lon;
            }
            taxonLocations[leftChild.getNr()][0] = taxonLocations[leftChild.getNr()][0] + lat;
            taxonLocations[leftChild.getNr()][1] = taxonLocations[leftChild.getNr()][1] + lon;
        }
        //find number of steps and changes of directions on the right branch
        numSteps = (int)(rTimeElapsed/timeStep);
        //find location of right child
        lat = 0.0;
        lon = 0.0;
        //set right child location initially to that of parent
        taxonLocations[rightChild.getNr()][0] = taxonLocations[root.getNr()][0];
        taxonLocations[rightChild.getNr()][1] = taxonLocations[root.getNr()][1];
        for (int i = 0; i < numSteps; i++){
            direction = Randomizer.nextDouble() * Math.PI;
            direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counterclockwise
            lat = Math.sin(direction)*spatialStep;
            lon = Math.cos(direction)*spatialStep;
            proposedLat = taxonLocations[rightChild.getNr()][0] + lat;
            proposedLon = taxonLocations[rightChild.getNr()][1] + lon;
            //boundary handling
            while (proposedLat > maxLatBoundary || proposedLat < minLatBoundary || proposedLon > maxLonBoundary || proposedLon < minLonBoundary){
                direction = Randomizer.nextDouble() *Math.PI;
                direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counter-clockwise
                lat = Math.sin(direction)*spatialStep;
                lon = Math.cos(direction)*spatialStep;
                proposedLat = taxonLocations[rightChild.getNr()][0] + lat;
                proposedLon = taxonLocations[rightChild.getNr()][1] + lon;
            }
            taxonLocations[rightChild.getNr()][0] = taxonLocations[rightChild.getNr()][0] + lat;
            taxonLocations[rightChild.getNr()][1] = taxonLocations[rightChild.getNr()][1] + lon;
        }
        if (leftChild.getChildCount() != 0){
            SRWWithBoundaries(leftChild);
        }
        if (rightChild.getChildCount() != 0){
            SRWWithBoundaries(rightChild);
        }
    }
}
