package bdtree.simulator;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.parser.XMLParserException;
import beast.base.util.Randomizer;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Description("Directly simulate trees from BDSS model.")
@Citation(value = "Stadler, T., & Yang, Z. (2013). Dating phylogenies with \n" +
        "  sequentially sampled tips. Systematic biology, 62(5), 674-688", DOI = "https://doi.org/10.1093/sysbio/syt030",
        year = 2012, firstAuthorSurname = "Stadler")

public class BirthDeathSerialSamplingTree extends beast.base.inference.Runnable {
    final public Input<TraitSet> traitSetInput = new Input<>("trait", "node times of sampled tips");
    final public Input<TaxonSet> taxonsetInput = new Input<>("taxonset", "set of taxa that correspond to the leafs in the tree");
    final public Input<Integer> speciesNrInput = new Input<>("speciesNr", "number of species in the simulated tree");
    final public Input<String> tipTimesInput = new Input<>("tipTime", "traits encoded as taxon=value pairs separated by commas");

    final public Input<String> m_outputFileNameInput = new Input<>("outputFileName", "If provided, simulated alignment is written to this file rather than to standard out.", Input.Validate.REQUIRED);
    final public Input<Integer> iterationsInput = new Input<>("iterations","number of trees to simulate", Input.Validate.REQUIRED);
    final public Input<Integer> logEveryInput = new Input<>("logEvery","frequency of screen log.", Input.Validate.REQUIRED);

    final public Input<RealParameter> rootAgeInput = new Input<>("rootAge", "the height of simulated trees", Input.Validate.REQUIRED);
    final public Input<RealParameter> birthRateInput = new Input<>("birthRate", "birth rate parameter", Input.Validate.REQUIRED);
    final public Input<RealParameter> deathRateInput = new Input<>("deathRate", "death rate parameter", Input.Validate.REQUIRED);
    final public Input<RealParameter> rhoInput = new Input<>("rho", "probability of sampling each extant lineage", Input.Validate.REQUIRED);
    final public Input<RealParameter> psiInput = new Input<>("psi", "sampling rate parameter, for analysis including fossils");


    //tree used for generating samples
    private Tree m_tree;
    private TraitSet traitSet;
    private TaxonSet taxonSet;
    private int speciesNr;
    private String traitName = "date-backward";

    // BDSS model parameters
    private double lambda;
    private double mu;
    private double rho;
    private double psi;
    private double tmrca;
    private int specifiedTreeNr;

    // name of output file
    private String m_outputFileName;
    private Integer logEvery = 10;

    // constants
    private double c1;
    private double c2;
    private double gt;

    private Node[] m_nodes;
    @Override
    public void initAndValidate() {
        // get inputs
        tmrca = rootAgeInput.get().getValue();

        if(taxonsetInput.get() == null){
            // if taxa names are not input, sp_i will be used
            speciesNr = speciesNrInput.get();
            taxonSet = new TaxonSet();
            taxonSet.initByName("taxon", creatTaxonSet());
        } else {
            taxonSet = taxonsetInput.get();
            speciesNr = taxonSet.getTaxonCount();
        }

        if(traitSetInput.get() == null){
            // if fossils are not input, random tip ages will be assigned
            traitSet = new TraitSet();
            StringBuilder tipAgesBuf =new StringBuilder();
            if(tipTimesInput.get() == null){
                for (int i = 0; i < (speciesNr - 1); i++){
                    tipAgesBuf.append(taxonSet.getTaxaNames().toArray()[i] + "=" + Randomizer.uniform(0, tmrca) + ",");
                }
                tipAgesBuf.append(taxonSet.getTaxaNames().toArray()[speciesNr-1] + "=" + Randomizer.uniform(0, tmrca));
                String randomTipAges = tipAgesBuf.toString();
                traitSet.initByName("value", randomTipAges , "taxa", taxonSet, "traitname", traitName);
            } else {
                traitSet.initByName("value", tipTimesInput.get(), "taxa", taxonSet, "traitname", traitName);
            }
        } else {
            traitSet = traitSetInput.get();
        }

        lambda = birthRateInput.get().getValue();
        mu = deathRateInput.get().getValue();
        rho = rhoInput.get().getValue();
        psi = psiInput.get().getValue();

        // check valid inputs
        if(lambda < mu) {
            throw new IllegalArgumentException("Birth rate is smaller than death rate.");
        }
        if (psi < 0 || (rho > 1 && rho < 0) || tmrca < 0) {
            throw new IllegalArgumentException("Illegal parameter values.");
        }

        // calculate the constants in the simulating functions
        c1 = Math.sqrt(Math.pow(lambda - mu - psi, 2) + 4 * lambda * psi);
        c2 = -(lambda - mu - 2 * lambda * rho - psi) / c1;
        gt = 1 / (Math.exp(-c1 * tmrca) * (1 - c2) + (1 + c2));

        // initialize to write the output .trees file
        m_tree = randomTreeTopology();
        m_tree.getRoot().setHeight(tmrca);

        m_outputFileName = m_outputFileNameInput.get();
        specifiedTreeNr = iterationsInput.get();
        logEvery = logEveryInput.get();

        Log.info("BDSSTreeSimulator:: Simulating trees from Birth-death process with sequentially sampled tips.");
    }

    @Override
    public void run() throws IllegalArgumentException, IllegalAccessException, IOException, XMLParserException {
        final long startTime = System.currentTimeMillis();

        // prepare the output .trees file
        PrintStream pstream;
        if (m_outputFileName == null) {
            pstream = System.out;
        }
        else {
            pstream = new PrintStream(m_outputFileName);
        }
        m_tree.init(pstream);
        pstream.print("\n");

        // simulating required number of trees
        for (int i = 0; i < specifiedTreeNr; i++) {
            // step1: initialize the tree
            m_tree = randomTreeTopology();
            m_tree.getRoot().setHeight(tmrca);

            // step2: simulating process
            drawDivTimes(m_tree);

            // step3: change tree for invalid node times
            m_nodes = m_tree.getNodesAsArray();
            // collect heights
            final double[] heights = new double[m_tree.getNodeCount()];
            final int[] reverseOrder = new int[m_tree.getNodeCount()];
            collectHeights(m_tree.getRoot(), heights, reverseOrder, 0);

            // adjust the tree topology to the divergence times
            final Node root = reconstructTree(heights, reverseOrder, 0, heights.length, new boolean[heights.length]);
            root.setParent(null);
            m_tree.setRoot(root);

            // step4: write simulated tree to the output file
            pstream.print("tree STATE_" + i + " = ");
            pstream.print(toNewick(m_tree.getRoot()));
            pstream.print(";" + "\n");

            // print screen log information
            if(i > 0 && i % logEvery == 0) {
                System.out.print(i + " trees successfully simulated." + "\n");
            }
        }

        // finish the .trees file
        m_tree.close(pstream);

        // screen info
        Log.info("BDSSTreeSimulator:: " + specifiedTreeNr + " were successfully simulated.");
        final long endTime = System.currentTimeMillis();
        Log.info.println("Total calculation time: " + (endTime - startTime) / 1000.0 + " seconds");
    }

    /*
     * Generate a random tree topology given the number of tips and tip ages
     */
    private Tree randomTreeTopology() {
        Tree ini_tree = new Tree();
        ini_tree.initByName("trait", traitSet, "taxonset", taxonSet);
        ini_tree.getRoot().setHeight(tmrca);

        // get the tip information to populate the random tree
        List<Node> species = ini_tree.getExternalNodes();
        List<Node> activeNodes = new ArrayList<>();
        int k = 0;
        for (Node node : species) {
            /*
             * Depend on Node.class in BEAST2 havin the following initiate method
             *     public Node(final String id, final double height) {
                    setID(id);
                    initAndValidate();
                    setHeight(height);
                   }
             */
            // no such constructor in beast 2.7
//            Node newNode = new Node(node.getID(), node.getHeight());
            Node newNode = new Node(node.getID());
            newNode.setHeight(node.getHeight());

            newNode.setNr(k);
            activeNodes.add(newNode);
            k += 1;
        }

        while (activeNodes.size() > 1) {

            // draw 2 random nodes from the list
            Node a = activeNodes.remove(Randomizer.nextInt(activeNodes.size()));
            Node b = activeNodes.remove(Randomizer.nextInt(activeNodes.size()));

            // add children to the parent
            Node parent = new Node();
            parent.addChild(a);
            parent.addChild(b);
            parent.setNr(k);

            k += 1;

            activeNodes.add(parent);
        }

        return new Tree(activeNodes.get(0));
    }

    /*
     * This method is called when taxonSet is not input
     * and therefore is used to creat a taxonSet by the specified number of species
     */
    private List<Taxon> creatTaxonSet(){
        Taxon taxon = new Taxon();
        List<String> taxaNames = new ArrayList<>(speciesNr);
        for (int i = 1; i < speciesNr + 1; i++){
            taxaNames.add("sp_" + i);
        }
        return taxon.createTaxonList(taxaNames);
    }

    private void drawDivTimes(Tree i_tree)  {
        int k;

        // iterate internal nodes except the root
        for (int j = speciesNr; j < i_tree.getRoot().getNr(); j++) {
            if (j != i_tree.getRoot().getNr()) {
                // step1: get z^* in Equation (4) in Stadler and Yang 2013
                // find tip on the left
                for (k = i_tree.getNode(j).getChild(1).getNr(); k >= i_tree.getLeafNodeCount(); k = i_tree.getNode(k).getChild(0).getNr());
                double z0 = i_tree.getNode(k).getHeight();
                // find tip on the right
                for (k = i_tree.getNode(j).getChild(0).getNr(); k >= i_tree.getLeafNodeCount(); k = i_tree.getNode(k).getChild(1).getNr());
                double z1 = i_tree.getNode(k).getHeight();
                double zstar = Math.max(z0, z1);

                // step2
                // calculate 1/g(z*)
                double gzstar = 1 / (Math.exp(-c1 * zstar) * (1 - c2) + (1 + c2));
                // a2 = (1/g(t_mrca)) - (1/g(z^*))
                double a2 = gt - gzstar;

                // step3
                // the constant part in the integral, which is H(zstar) and H is CDF of divergence times
                double constantChildren = 1.0 / (a2 * (((1 - c2) * Math.exp(-c1 * zstar)) + (1 + c2)));

                // step4: drawn a random number of Uniform(0,1)
                double y = Randomizer.uniform(0, 1);

                // step5: calculate the inverse function, i.e. H^(-1)
                double x = Math.log((1 / (a2 * (y + constantChildren) * (1 - c2))) - ((1 + c2) / (1 - c2))) / (-c1);

                // step6: set the simulated divergence time
                i_tree.getNode(j).setHeight(x);
            }
        }
    }


    // write BEAST tree to newick format in the output .trees file
    private String toNewick(Node node) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft()));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight()));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr()+1);
        }

        buf.append(":");
        appendDouble(buf, node.getLength());
        return buf.toString();
    }


    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        buf.append(d);
    }


    private Node reconstructTree(final double[] heights, final int[] reverseOrder, final int from, final int to, final boolean[] hasParent) {
        //nodeIndex = maxIndex(heights, 0, heights.length);
        int nodeIndex = -1;
        double max = Double.NEGATIVE_INFINITY;
        for (int j = from; j < to; j++) {
            if (max < heights[j] && !m_nodes[reverseOrder[j]].isLeaf()) {
                max = heights[j];
                nodeIndex = j;
            }
        }
        if (nodeIndex < 0) {
            return null;
        }
        final Node node = m_nodes[reverseOrder[nodeIndex]];

        //int left = maxIndex(heights, 0, nodeIndex);
        int left = -1;
        max = Double.NEGATIVE_INFINITY;
        for (int j = from; j < nodeIndex; j++) {
            if (max < heights[j] && !hasParent[j]) {
                max = heights[j];
                left = j;
            }
        }

        //int right = maxIndex(heights, nodeIndex+1, heights.length);
        int right = -1;
        max = Double.NEGATIVE_INFINITY;
        for (int j = nodeIndex + 1; j < to; j++) {
            if (max < heights[j] && !hasParent[j]) {
                max = heights[j];
                right = j;
            }
        }

        node.setLeft(m_nodes[reverseOrder[left]]);
        node.getLeft().setParent(node);
        node.setRight(m_nodes[reverseOrder[right]]);
        node.getRight().setParent(node);
        if (node.getLeft().isLeaf()) {
            heights[left] = Double.NEGATIVE_INFINITY;
        }
        if (node.getRight().isLeaf()) {
            heights[right] = Double.NEGATIVE_INFINITY;
        }
        hasParent[left] = true;
        hasParent[right] = true;
        heights[nodeIndex] = Double.NEGATIVE_INFINITY;


        reconstructTree(heights, reverseOrder, from, nodeIndex, hasParent);
        reconstructTree(heights, reverseOrder, nodeIndex, to, hasParent);
        return node;
    }

    private int collectHeights(final Node node, final double[] heights, final int[] reverseOrder, int current) {
        if (node.isLeaf()) {
            heights[current] = node.getHeight();
            reverseOrder[current] = node.getNr();
            current++;
        } else {
            current = collectHeights(node.getLeft(), heights, reverseOrder, current);
            heights[current] = node.getHeight();
            reverseOrder[current] = node.getNr();
            current++;
            current = collectHeights(node.getRight(), heights, reverseOrder, current);
        }
        return current;
    }
}
