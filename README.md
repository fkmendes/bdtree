# bdtree

A central repository for birth-death models in BEAST 2.

**Contributors**:   Fábio K. Mendes, Rong Zhang

## Building

In order to build the *bdtree* package .jar file:

(1) Clone and build BEAST2 from [here](https://github.com/CompEvol/beast2);

(2) Clone the *bdtree* repository side-by-side with the *beast2/* directory resulting from step (1);

(3) From *bdtree/*, type:

```
$ ant
```

This command invokes the *ant* tool, which executes the building instructions inside *build.xml*. The new *bdtree* .jar file will be put inside *build/dist/*.

## Running an analysis with *bdtree*

From your shell or Terminal, type:

```
$ /path/to/build/dist/bdtree.v0.0.1.jar /path/to/an_analysis.xml
```

## Models

Below we list all the models contained within *bdtree*:

**Birth-death-sequential-sampling (BDSS)**

[[ref]](https://academic.oup.com/sysbio/article/62/5/674/1684217) Tanja Stadler, Ziheng Yang (2013). Dating phylogenies with sequentially sampled tips. *Syst. Biol.* 62(5), 674-688.
