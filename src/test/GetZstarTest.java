package test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.util.TreeParser;
import org.junit.Assert;
import org.junit.Test;
import bdtree.likelihood.BirthDeathSequentialSampling;
import beast.evolution.tree.Node;
import java.util.ArrayList;
import java.util.List;

public class GetZstarTest {

    @Test
    public void testingGetZStar() {
        // (1) set up the test
        // define a tree
        String treeStr = "((A:2.0,B:6.0):1,((C:2.0,D:1.0):3.0,(E:1.0,(F:2.0,G:2.0):2.0):1.0):2.0):0.0;";

        List<Taxon> taxaSet = new ArrayList<>();
        taxaSet.add(0, new Taxon("A"));
        taxaSet.add(1, new Taxon("B"));
        taxaSet.add(2, new Taxon("C"));
        taxaSet.add(3, new Taxon("D"));
        taxaSet.add(4, new Taxon("E"));
        taxaSet.add(5, new Taxon("F"));
        taxaSet.add(6, new Taxon("G"));
        TaxonSet taxonSet = new TaxonSet();
        taxonSet.initByName("taxon", taxaSet);

        // specify sampled dates
        String fossilAges = "A=4.0,B=0.0,C=0.0,D=1.0,E=3.0,F=0.0,G=0.0";
        String TraitName = "date-backward";
        TraitSet fossilSet = new TraitSet();
        fossilSet.initByName("value", fossilAges,"traitname", TraitName, "taxa", taxonSet);

        TreeParser tree = new TreeParser();
        tree.initByName("newick", treeStr, "taxonset", taxonSet, "trait", fossilAges, "IsLabelledNewick", true, "adjustTipHeights", false);

        // define BDSS model
        BirthDeathSequentialSampling treePrior = new BirthDeathSequentialSampling();
        RealParameter birthRate = new RealParameter(new Double[] {1.0});
        RealParameter deathRate = new RealParameter(new Double[] {1.0});
        RealParameter samplingRate = new RealParameter(new Double[] {0.001});
        RealParameter rho = new RealParameter(new Double[] {0.0});
        treePrior.initByName("tree", tree, "birthRate", birthRate, "deathRate", deathRate, "psi", samplingRate, "rho", rho, "rootAge", 7.0);

        // get the internal nodes of the tree
        List<Node> internalNodes = tree.getInternalNodes();

        // (2) run the test
        // zstar for (A, B) = max(4, 0)
        Assert.assertEquals(4.0, treePrior.getzstar(internalNodes.get(0)), 1e-4);

        // zstar for (C, D) = max(0, 1)
        Assert.assertEquals(1.0, treePrior.getzstar(internalNodes.get(1)), 1e-4);

        // zstar for (G, F) = max(0, 0)
        Assert.assertEquals(0.0, treePrior.getzstar(internalNodes.get(2)), 1e-4);

        // zstar for (F, E) = max(0, 3)
        Assert.assertEquals(3.0, treePrior.getzstar(internalNodes.get(3)), 1e-4);

        // zstar for (D, E) = max(1, 3)
        Assert.assertEquals(3.0, treePrior.getzstar(internalNodes.get(4)), 1e-4);
    }
}