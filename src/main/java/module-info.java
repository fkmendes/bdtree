open module bdtree {
    requires beast.base;
    requires beast.pkgmgmt;

    exports bdtree.likelihood;
    exports bdtree.simulator;

    provides beast.base.core.BEASTInterface with
        bdtree.likelihood.BirthDeathSequentialSampling,
        bdtree.simulator.BirthDeathSerialSamplingTree;
}
