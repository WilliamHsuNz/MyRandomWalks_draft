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
import beast.util.Randomizer;
import beast.util.XMLProducer;
import sphericalGeo.TreeTraitMap;

/**
 * Created by williamhsu on 18/03/16.
 */
@Description("Traits containing location data generated using " +
        "simple isotropic random walk are mapped onto a given tree.")
public class SimulatedTreeTraitMap extends TreeTraitMap {
    //final public Input<TreeTraitMap> m_dataInput = new Input<>("data", "trait data which specifies... ", Validate.REQUIRED);
    //final public Input<Tree>m_treeInput = new Input<Tree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    TreeInterface tree;
    RealParameter rP;
    String traitName;
    double [][] taxonLocations;
    double timeStep;
    double spatialStep;
    double direction;
    //TreeTraitMap ttm;

    public SimulatedTreeTraitMap(){
        parameterInput.setRule(Validate.OPTIONAL);
    }

    public void initAndValidate(){
        //ttm = m_dataInput.get();
        //ttm = m_data.get();
        tree = treeInput.get();
        rP = parameterInput.get();

        if (!(tree instanceof Tree))
            throw new IllegalArgumentException("Tree input must be a true Tree, not just TreeInterface.");

        makeRootLocation((Tree)tree);//this needs to be a Tree object it is currently a tree interface

       // traitName = traitNameInput.get();
        super.initAndValidate();
       // System.out.println(ttm.getTraitName());
        //getLeafNodeCount());
    }

    private void makeRootLocation(Tree tree){
        Node root = tree.getRoot();
        System.out.println("Tree Height: " + root.getHeight());

        double [] rootLocation = {0.0, 0.0};

        taxonLocations = new double [tree.getNodeCount()][2];
        taxonLocations[root.getNr()] = rootLocation;

        //double angle = 0.0;
        SRW(root);
    }

    private void SRW(Node root){
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
            direction = Randomizer.nextDouble() * Math.PI;
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
            direction = Randomizer.nextDouble() * Math.PI;
            direction = (Randomizer.nextBoolean() == true)? direction : direction* -1;//clockwise or counterclockwise
            lat = Math.sin(direction)*spatialStep;
            lon = Math.cos(direction)*spatialStep;
            taxonLocations[rightChild.getNr()][0] = taxonLocations[rightChild.getNr()][0] + lat;
            taxonLocations[rightChild.getNr()][1] = taxonLocations[rightChild.getNr()][1] + lon;
        }
        if (leftChild.getChildCount() != 0){
            SRW(leftChild);
        }
        if (rightChild.getChildCount() != 0){
            SRW(rightChild);
        }
    }

}