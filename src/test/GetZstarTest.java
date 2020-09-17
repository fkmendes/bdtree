package test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import org.junit.Assert;
import org.junit.Test;
import bdtree.likelihood.BirthDeathSequentialSampling;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import java.util.ArrayList;
import java.util.List;

public class GetZstarTest {

    /*
     * This test checks the value of Z star
     * for all internal nodes in a non-ultrametric
     * tree without sampled ancestors (those are
     * not supported by this likelihood)
     */
    @Test
    public void testingGetZStar() {

        String treeStr = "((A:2.0,B:6.0):1,((C:2.0,D:1.0):3.0,(E:1.0,(F:2.0,G:2.0):2.0):1.0):2.0):0.0;";

        List<Taxon> taxaSet = new ArrayList<>();
        taxaSet.add(new Taxon("A"));
        taxaSet.add(new Taxon("B"));
        taxaSet.add(new Taxon("C"));
        taxaSet.add(new Taxon("D"));
        taxaSet.add(new Taxon("E"));
        taxaSet.add(new Taxon("F"));
        taxaSet.add(new Taxon("G"));
        TaxonSet taxonSet = new TaxonSet();
        taxonSet.initByName("taxon", taxaSet);

        // specify sampled dates
        String fossilAges = "A=4.0,B=0.0,C=0.0,D=1.0,E=3.0,F=0.0,G=0.0";
        TraitSet fossilSet = new TraitSet();
        fossilSet.initByName("value", fossilAges,
                "traitname", "date-backward",
                "taxa", taxonSet);

        // initializing tree
        Tree tree = new TreeParser();
        tree.initByName("newick", treeStr,
                "taxonset", taxonSet,
                "trait", fossilAges,
                "IsLabelledNewick", true,
                "adjustTipHeights", false);

        // initializing BDSS likelihood
        BirthDeathSequentialSampling treePrior = new BirthDeathSequentialSampling();
        RealParameter birthRate = new RealParameter(new Double[] { 1.0 });
        RealParameter deathRate = new RealParameter(new Double[] { 1.0 });
        RealParameter samplingRate = new RealParameter(new Double[] { 0.001 });
        RealParameter rho = new RealParameter(new Double[] { 0.0 });
        treePrior.initByName("tree", tree,
                "birthRate", birthRate,
                "deathRate", deathRate,
                "psi", samplingRate,
                "rho", rho,
                "rootAge", 7.0);

        // grabbing all zStars from internal nodes
        Double[] zStars = new Double[tree.getInternalNodeCount() - 1]; // -1 b/c ignoring root
        int j = 0;
        for (Node i: tree.getInternalNodes()) {
            if (!i.isRoot()) {
                zStars[j] = treePrior.getZStar(i);
                j++;
            }
        }

        // test!
        Assert.assertArrayEquals(new Double[] { 4.0, 1.0, 0.0, 3.0, 3.0 }, zStars);
    }
}