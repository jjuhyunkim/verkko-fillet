## Preprocessing : `pp`

```{eval-rst}
.. module:: verkkofillet.pp
```

```{eval-rst}
.. currentmodule:: verkkofillet
```

Extracting gap from path, preprocessing recipes. 

### Reading and writing
```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.read_Verkko
   pp.readChr
   pp.writeFixedPaths
   pp.save_Verkko
   pp.load_Verkko
   pp.mkCNSdir
   pp.updateCNSdir_missingEdges
   pp.loadGiraffe
   pp.readGaf
   pp.getQV

``` 

### Preprocessing
```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.find_intra_telo
   pp.find_reads_intra_telo
   pp.findGaps
   pp.searchSplit
   pp.searchNodes
   pp.highlight_nodes
   pp.fillGaps
   pp.estLoops
   pp.checkGapFilling
   pp.naming_contigs
   pp.find_multi_used_node
   pp.find_hic_support

```

### Last touch
```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: generated/

   pp.get_NodeChr
   pp.find_multiContig_chr

   
```