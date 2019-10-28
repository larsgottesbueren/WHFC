# WHFC
Implementation of the Weighted HyperFlowCutter hypergraph partitioning algorithm.
It works as a stand-alone partitioner but it works best as a refinement algorithm on a given partition.
This is a header-only library, so no setup overhead, if you want to use it in your own tool.
If you want an example integration, check out the following links [flow hypergraph extraction](https://github.com/larsgottesbueren/kahypar/blob/HyperFlowCutterIntegration/kahypar/partition/refinement/flow/whfc_flow_hypergraph_extraction.h), [2-way refinement](https://github.com/larsgottesbueren/kahypar/blob/HyperFlowCutterIntegration/kahypar/partition/refinement/flow/2way_hyperflowcutter_refiner.h) and [k-way refinement](https://github.com/larsgottesbueren/kahypar/blob/HyperFlowCutterIntegration/kahypar/partition/refinement/flow/kway_hyperflowcutter_refiner.h) for my integration in [KaHyPar](https://github.com/SebastianSchlag/kahypar/tree/master/kahypar). 
