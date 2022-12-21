package bdtree.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.Distribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;

import java.util.List;
import java.util.Random;

@Description("This class implements BDSS model for tree likelihood that is used in mcmcTree. " +
        "If only extant species are analyzed, the class returns the density of the birth-death process " +
        "described in Yang and Rannala 2006. If fossil species are included, the class returns the density of " +
        "birth-death with sequential sampling process described in Stader and Yang 2013.")
// Note that density for the calibration of divergence times (internal nodes) has not been implemented yet.
// Note: to comment tree.getRoot.sort() in tree logger
public class BirthDeathSequentialSampling extends Distribution {
    final public Input<Tree> treeInput = new Input<>("tree", "tree over which to calculate a prior or likelihood");
    final public Input<RealParameter> birthRateInput = new Input<>("birthRate", "birth rate parameter", Input.Validate.REQUIRED);
    final public Input<RealParameter> deathRateInput = new Input<>("deathRate", "death rate parameter", Input.Validate.REQUIRED);
    final public Input<RealParameter> rhoInput = new Input<>("rho", "probability of sampling each extant lineage", Input.Validate.REQUIRED);
    final public Input<RealParameter> psiInput = new Input<>("psi", "sampling rate parameter, for analysis including fossils");
    final public Input<Double> rootAgeLowerInput = new Input<>("lower", "lower soft bound for the root age");
    final public Input<Double> rootAgeUpperInput = new Input<>("upper", "upper soft bound for the root age");
    final public Input<Double> rootAgeInput = new Input<>("rootAge", "specified root age when the tree height is fixed");

    private double birthRate;
    private double deathRate;
    private double rho;
    private double psi;
    private double lower;
    private double upper;

    private double tmrca;
    private Tree tree;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        // get input parameters
        getBDSSModelParameters();
        if (birthRate < deathRate) {
            throw new IllegalArgumentException("Birth rate is smaller than death rate.");
        }
        if (psi < 0) {
            throw new IllegalArgumentException("Sampling rate is negative.");
        }

        // get sampled tree
        tree = treeInput.get();

        // get root age or t_mrca
        if (rootAgeInput.get() == null) {
            // if rootAge is not input, the tree height is sampled under a Uniform prior with soft boundaries
            lower = rootAgeLowerInput.get();
            upper = rootAgeUpperInput.get();
            if (lower >= upper) {
                throw new IllegalArgumentException("The bounds for root age is not correctly specified.");
            }
        } else {
            // if the rootAge is input, the tree height will be fixed
            tmrca = rootAgeInput.get();
            tree.getRoot().setHeight(tmrca);
        }
    }

    private void getBDSSModelParameters(){
        birthRate = birthRateInput.get().getValue();
        deathRate = deathRateInput.get().getValue();
        rho = rhoInput.get().getValue();
        if (psiInput.get() == null) {
            psi = 0.0;
        } else {
            psi = psiInput.get().getValue();
        }
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        // get BDSS model parameters
        getBDSSModelParameters();

        // get the root age, i.e. height of the tree
        double t1 = tmrca;
        if (rootAgeInput.get() == null) {
            t1 = tree.getRoot().getHeight();
            // Add the density of the root age
            logP += lnptCalibrationDensity(t1, lower, upper);
        } else {
            // The tree height is fixed
            tree.getRoot().setHeight(tmrca);
        }

        List<Node> internalNodes = tree.getInternalNodes();
        // This calculates the prior density of node times expect the root in the species tree.
        logP += lnpriorTimesBDS_Approach1(internalNodes, t1, birthRate, deathRate, rho, psi);

        return logP;
    }

    /*
     * This method implements the same method from mcmctree.c
     * which gives that probability density of divergence times given the sampling times and root age.
     *
     */
    private double lnpriorTimesBDS_Approach1(List<Node> internalNodes, double t1, double lambda, double mu, double rho, double psi) {
        double logP = 0.0;

        /* psi = 0  Birth-death with species sampling (Yang & Rannalar 2006)  */
        // Case 1: birth rate is equal to death rate, no fossils included.
        // This case is explained by Equation (8) in Yang & Rannalar 2006, i.e.
        // g(t) = (1 + rho * lambda * t1) / [t1 * (1 + rho * lambda * t)^2]
        // Both numerator and denominator are divided by t1,
        // we get numerator denoted by c1 and denominator denoted by c2 * c2.
        if (lambda == mu && psi == 0.0) {
            double c1 = 1 / t1 + rho * lambda;
            // iterating all internal nodes expect root
            for (Node intNode : internalNodes)
                if (!intNode.isRoot()) {
                    double c2 = 1 + rho * lambda * intNode.getHeight();
                    logP += Math.log(c1 / (c2 * c2));
                }
        }

        // Case 2: birth rate is NOT equal to death rate, no fossils included.
        // This is the special case explained by Equation (6) in Stadler & Yang (2013 Syst Biol).
        else if (lambda != mu && psi == 0.0) {
            // (1) calculating the second product in Equation (6),
            //     which is constant over the following loop
            // a = [lambda * (1 - rho) - mu]
            double a = lambda - rho * lambda - mu;

            // e = exp[(mu - lambda) * t1]
            double e = Math.exp((mu - lambda) * t1);

            // the constant part
            double c1 = (rho * lambda + a * e) / (1 - e);

            if (Math.abs(1 - e) < 1E-100) {
                System.out.print("e too close to 1..");
            }

            // (2) iterating all internal nodes expect root
            for (Node intNode : internalNodes)
                if (!intNode.isRoot()) {
                    // e = exp[(mu - lambda) * t]
                    e = Math.exp((mu - lambda) * intNode.getHeight());

                    // c2 is the part that has a square, i.e.
                    // c2 = (lambda - miu) / {rho * lambda + [lambda - lambda * rho) - miu] * exp[(mu - lambda) * t}
                    double c2 = (lambda - mu) / (rho * lambda + a * e);

                    // here c2 combines the first and second product in Equation (6), i.e.
                    // c2^2 * exp[(mu - lambda) * t] * c1
                    c2 *= c2 * e * c1;

                    if (c2 < 1E-300 || c2 > 1E300) {
                        //System.out.print("c2 not numeric.");
                        return Double.NEGATIVE_INFINITY;
                    }

                    // accumulate the log density
                    logP += Math.log(c2);
                }
        }

        /* psi != 0   Birth-death with sequential sampling model (Stadler & Yang 2013) */
        // Case 3: for analysis including fossils
        // The calculation is detailed by Equation (4) in Stadler & Yang (2013 Syst Biol).
        else {
            // c1 = square root of (lambda - mu - psi)^2 + 4 * lambda * psi
            double c1 = Math.sqrt(Math.pow(lambda - mu - psi, 2) + 4 * lambda * psi);

            double c2 = -(lambda - mu - 2 * lambda * rho - psi) / c1;

            // g(t) = [exp(-c1 * t) * (1 - c2)] + (1 + c2)
            // gt1 calculates 1/g(t1).
            double gt1 = 1 / (Math.exp(-c1 * t1) * (1 - c2) + (1 + c2));

            if (gt1 < -1E-300 || gt1 > 1E300) {
                System.out.print("gt1 not numeric.");
            }

            // iterating all internal nodes expect root
            // here we iterate by node index, i.e. excluding the index of tip
            // iterate internal nodes except the root
            for (Node aNode : internalNodes) {
                if (!aNode.isRoot()) {
                    // get the older sampled tips below this internal node
                    double zstar = getZStar(aNode);

                    // calculate 1/g(z*)
                    zstar = 1 / (Math.exp(-c1 * zstar) * (1 - c2) + (1 + c2));

                    double t = aNode.getHeight();
                    // calculate g(t) for this internal node
                    double gt = Math.exp(-c1 * t) * (1 - c2) + (1 + c2);

                    if(gt1 == zstar) {
                        return Double.NEGATIVE_INFINITY;
                    }

                    // accumulate the log density
                    // numerator = c1 * (1 - c2 ) * exp(-c1 * t)
                    // denominator = g(t)^2 * [1/g(t1) - 1/g(z*)]
                    logP += -c1 * t + Math.log(c1 * (1 - c2) / (gt * gt * (gt1 - zstar)));
                }
            }
        }

        return logP;
    }

    /*
     * The method calculates the log density of root age with soft bound.
     * The calculation is detailed in Equation (17) in Yang and Rannala 2006.
     * t: root age
     * tL: lower bound
     * tU: upper bound
     */
    private double lnptCalibrationDensity(double t, double tL, double tU) {
        // probability that the root age exceeds the lower and upper bounds
        // currently, the default value is 0.025
        double tailL = 0.025;
        double tailU = 0.025;

        // case 1: the root age is within the the range
        if (t > tL && t < tU) {

            return Math.log((1 - tailL - tailU) / (tU - tL));
        }

        // case 2: the root age is smaller than lower bound
        else if (t < tL) {
            // theta1 = 0.95 * tL / [0.025 * (tL - tR)]
            double theta1 = (1 - tailL - tailU) * tL / (tailL * (tU - tL));

            // f(t) = [0.025 * (theta1 / tL)] * (t / tL)^(theta1 - 1)
            // log(f(t)) = log( 0.025 * theta1 / tL) + (theta1 - 1) * log(t / tL)
            return Math.log(tailL * theta1 / tL) + (theta1 - 1) * Math.log(t / tL);
        }

        // case 3: the root age is larger than upper bound
        else {
            // theta2 = 0.95 / [0.025 * (tL - tR)]
            double theta2 = (1 - tailL - tailU) / (tailU * (tU - tL));

            // f(t) = (0.025 * theta2) * exp[-theta2 * (t - tU)]
            // log(f(t)) = log(0.025 * theta2) + [-theta2 * (t - tU)]
            return Math.log(tailU * theta2) - theta2 * (t - tU);
        }
    }

    // for Junit test
    public double getZStar(Node aNode) {
        Node child;
        // step1: get z^* in Equation (4) in Stadler and Yang 2013
        // find tip on the left
        for (child = aNode.getChild(1); !child.isLeaf(); child = child.getChild(0)) ;
        double z0 = child.getHeight();
        // find tip on the right
        for (child = aNode.getChild(0); !child.isLeaf(); child = child.getChild(1)) ;
        double z1 = child.getHeight();


        // compare the node times of the two tips
        return Math.max(z0, z1);
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() || birthRateInput.get().somethingIsDirty() ||
                InputUtil.isDirty(deathRateInput) || InputUtil.isDirty(psiInput) ||
                InputUtil.isDirty(rhoInput);
    }

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub

    }
}

