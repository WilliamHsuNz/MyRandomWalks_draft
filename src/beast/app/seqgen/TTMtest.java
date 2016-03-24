package beast.app.seqgen;

import beast.core.BEASTInterface;
import beast.core.parameter.RealParameter;
import ec.util.MersenneTwisterFast;
import sim.util.distribution.VonMises;
import sphericalGeo.TreeTraitMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.List;
import java.util.Arrays;

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
 * Created by williamhsu on 17/03/16.
 */
public class TTMtest {
    double [] a = new double[2];
    //a[0] = 0.4;
    //a[1] = 0.3;

    //RealParameter rP = new RealParameter(a);



    public static void main(String[] args) {
        TTMtest test = new TTMtest();
    }
}
